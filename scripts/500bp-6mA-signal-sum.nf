#!/usr/bin/env/ nextflow
nextflow.enable.dsl = 2
pod5_dir = "pod5_pass/"
dorado_dir = "/user/home/app/dorado-version-linux-x64/bin"
model = "/path/to/model"
reference = "path/to/reference.fa"
enzyme = "MboI"
prefix = "sample-name"
process basecalling {
    input:
        tuple val(sample_id), path(pod5_dir)
    output:
        path("${sample_id}.bam")
    """
    ${dorado_dir}/dorado basecaller ${dorado_basecall_model} ${pod5_dir} --modified-bases 6mA 5mC_5hmC | samtools view -bhS - > ${sample_id}.bam
    """
}
process porec_digest {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        path("${sample_id}.${enzyme}.digest.bam")
    script:
    """
    pore-c-py digest ${bam_file} ${enzyme} --output ${sample_id}.${enzyme}.digest.bam
    """
}
process mapping {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        path("${sample_id}.dorado.map.bam")
    script:
    """
    ${dorado_dir}/dorado aligner ${reference} ${sample_id}.${enzyme}.digest.bam > ${sample_id}.dorado.map.bam
    """
}
process first-filtering {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        path("${sample_id}.dorado.filter.bam")
    script:
    """
    samtools view -b -h -F 4 ${sample_id}.dorado.map.bam >${sample_id}.dorado.filter.bam
    """
}
process second-filtering {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        path("${sample_id}.methylation.bam")
    script:
    """
    samtools view -hb -q 30 ${sample_id}.dorado.filter.sort.bam -o ${sample_id}.methylation.bam
    """
}
process sorting {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        path("${sample_id}.methylation.sort.bam")
    script:
    """
    samtools sort ${sample_id}.methylation.bam >${sample_id}.methylation.sort.bam
    """
}
process index {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        tuple val(sample_id), path("${sample_id}.methylation.sort.bam.bai")
    script:
    """
    samtools index ${sample_id}.methylation.sort.bam
    """
}
process extract_modification {
    input:
        tuple val(sample_id), path(bam_file)
    output:
        tuple val(sample_id), path("${sample_id}.6mA.pass.bedgraph")
    script:
    """
    modbam2bed -m 6mA ${reference.fa} ${bam_file} > ${sample_id}.6mA.bed
    awk '{print \$1,\$2,\$3,\$11}' ${sample_id}.6mA.bed | awk '\$4 > 0 && \$4 != "nan"' > ${sample_id}.6mA.pass.bedgraph
    """
}
process index_fasta {
    input:
    path(fasta)
    output:
    path("${reference.fa}.fai")
    script:
    """
    samtools faidx ${reference.fa}
    """
}
process extract_chrom_size {
    input:
    path(fai_file)
    output:
    path("reference.chrom.size")
    script:
    """
    cut -f1,2 ${reference.fa}.fai > reference.chrom.size
    """
}
process make_windows {
    input:
    path(chrom_size_file)
    output:
    path("reference.chrom.size.500")
    script:
    """
    bedtools makewindows -g reference.chrom.size -w 500 > reference.chrom.size.500
    """
}
process sum_6mA_signal {
    input:
    path(window_bed)
    path(bedgraph_file)
    output:
    path("${params.prefix}.500bp.6mA.sum.txt")
    script:
    """
    bedtools intersect -a reference.chrom.size.500 -b ${sample_id}.6mA.pass.bedgraph -wa -wb >500bp.6mA.txt
    awk '{key=\$1"\\t"\$2"\\t"\$3; sum[key]+= \$7} END {for (k in sum) print k, sum[k]}' 500bp.6mA.txt | sort -k1,1 -k2,2n -k3,3n > 500bp.6mA.sum.txt
    """
}
workflow {
    pod5_samples | basecalling \
                 | porec_digest \
                 | mapping \
                 | first_filtering \
                 | sorting \
                 | index \
                 | extract_modification \
                 | sum_6mA_signal(make_windows(extract_chrom_size(genome_index(params.reference).out).out).out)
}
