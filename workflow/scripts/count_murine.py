import csv  
import os                                                                                                                                                
                                                                                                                                                             
sample = snakemake.wildcards.sampleid                                                                                                                          
input_files = snakemake.input
                          
# Create output directories if they don't exist
for f in snakemake.output:
    os.makedirs(os.path.dirname(f), exist_ok=True)    

human_reads, mouse_reads, amb_reads = [], [], []

for file_path in input_files:
  with open(file_path, 'r') as file:
      next(file)  # skip header
      for line in file:
          if not line.strip():
              continue
          parts = line.strip().split('\t')
          header = parts[0]
          human_count = int(parts[1])
          mouse_count = int(parts[2])

          row = [header, human_count, mouse_count]                                                                                                           

          if human_count >= 1:                                                                                                                               
              human_reads.append(row)
          elif mouse_count >= 1:
              mouse_reads.append(row)
          else:
              amb_reads.append(row)
                                                                                                                                                             
# Count how many reads in each category
category_counts = {                                                                                                                                            
  "Human": len(human_reads),
  "Mouse": len(mouse_reads),                                                                                                                                 
  "Ambiguous": len(amb_reads),
}                                                                                                                                                              
              
# Write counts summary
with open(snakemake.output[0], 'w', newline='') as f:
  writer = csv.writer(f, delimiter='\t')
  writer.writerow(["Header", "Human", "Mouse", "Ambiguous"])
  writer.writerow([sample] + list(category_counts.values()))                                                                                                 

# Write detailed TSVs for each category                                                                                                                        
header = ["Header", "Human", "Mouse"]
category_data = {                                                                                                                                              
  snakemake.output[1]: human_reads,                                                                                                                          
  snakemake.output[2]: mouse_reads,
  snakemake.output[3]: amb_reads,                                                                                                                            
}               
                                                                                                                                                             
for output_file, data in category_data.items():
  with open(output_file, 'w', newline='') as f:
      writer = csv.writer(f, delimiter='\t')
      writer.writerow(header)                                                                                                                                
      writer.writerows(data)