#!/usr/bin/env/ nextflow
nextflow.enable.dsl = 2
params.cutoff = the-cutoff-value-calculated-by-R
params.sample_id = "sample-name" 
params.input_bam = "${params.sample_id}.methylation.bam"
params.sum_file = "500bp.6mA.sum.txt"
params.combine_script = "combine-6mA-get.py"
process get6mA {
    input:
    path params.sum_file
    output:
    path "6mA-get-sequence.txt"
    script:
    """
    awk '$4>${params.cutoff}' ${params.sum_file} > 6mA-get-sequence.txt
    """
}
process bam2sam {
    input:
    path params.input_bam
    output:
    path "${params.sample_id}.methylation.sam"
    script:
    """
    samtools view ${params.input_bam} | cut -f 1,3,4 > ${params.sample_id}.methylation.sam
    """
}
process sam2bed {
    input:
    path "${params.sample_id}.methylation.sam"
    output:
    path "${params.sample_id}.fragment.bed"
    script:
    """
    sed 's/:/\t/g' ${params.sample_id}.methylation.sam | \
    awk '{print $4, $5, $5+$3-$2, $1}' OFS='\\t' > ${params.sample_id}.fragment.bed
    """
}
process intersect6mA {
    input:
    path "6mA-get-sequence.txt"
    path "${params.sample_id}.fragment.bed"
    output:
    path "6mA-overlap-frag.txt"
    script:
    """
    bedtools intersect -a 6mA-get-sequence.txt -b ${params.sample_id}.fragment.bed -wa -wb | \
    awk '{print $8, $1, $2, $3}' | sort | uniq > 6mA-overlap-frag.txt
    """
}
process combine6mA {
    input:
    path "6mA-overlap-frag.txt"
    path params.combine_script
    output:
    path "6mA-interact.txt"
    script:
    """
    python ${params.combine_script}
    """
}
process filterInteraction {
    input:
    path "6mA-interact.txt"
    output:
    path "6mA-interact.10k.bedpe"
    script:
    """
    sort 6mA-interact.txt | uniq | awk '$5-$3>=10000' > 6mA-interact.10k.bedpe
    """
}
workflow {
    get6mA()
    bam2sam()
    sam2bed(bam2sam.out)
    intersect6mA(get6mA.out, sam2bed.out)
    combine6mA(intersect6mA.out)
    filterInteraction(combine6mA.out)
}
