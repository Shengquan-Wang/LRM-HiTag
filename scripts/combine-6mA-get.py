def generate_combinations(file_path, output_path):   
    group_data = {}  
    with open(file_path, 'r') as file:  
        for line in file:  
            columns = line.strip().split()  
            if len(columns) != 4:  
                continue    
            group = columns[0]   
            region = tuple(columns[1:4])  
            if group not in group_data:  
                group_data[group] = []  
            group_data[group].append(region)    
    results = []  
    for group, regions in group_data.items():  
        num_regions = len(regions)  
        for i in range(num_regions):  
            for j in range(i + 1, num_regions):  
                results.append(f"{regions[i][0]}\t{regions[i][1]}\t{regions[i][2]}\t{regions[j][0]}\t{regions[j][1]}\t{regions[j][2]}")  
    with open(output_path, 'w') as output_file:  
        output_file.write('\n'.join(results) + '\n')  
    print(f" The result has been stored {output_path}")  
input_file_path = '6mA-overlap-frag.txt'    
output_file_path = '6mA-interact.txt'  
generate_combinations(input_file_path, output_file_path)
