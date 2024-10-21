import itertools 
import random 
import pandas as pd

### references needed for pipeline
# hg38 = "/insert/path/to/Homo_sapiens_assembly38.fasta"
# vcf = "/insert/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"
hg38 = "/tgen_labs/barthel/references/GRCh38/hg38_all_viral.fasta"
vcf = "/tgen_labs/barthel/references/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf"
workdir: "/scratch/smankame/"  ##this is where all the slurm files will be created

### this path is where all your files will be made
base_path = "/scratch/smankame/"

### user makes a csv file where the first column is study name, the second column is sample name and the third and fourth columns are fastq_R1 and fastq_R2 files.
### study name (first column) will be used to name the folders and sample name (second column) will be used to name all subsequent files
sWGS_table = pd.read_csv("/tgen_labs/barthel/software/github/barthel/cfDNA/sWGS/config/sWGS_data_table.txt",sep='\t', header=None, names=["Study", "Sample", "R1","R2"])
sWGS_table.index = sWGS_table['Sample']

### if there are multiple studies within the csv table, study filter will work with specific subsets of samples. Change the study filter value to the study name. 
study_filter = "BCOH_UPN210"
filtered_sWGS_table = sWGS_table[sWGS_table['Study'] == study_filter]
study_id =pd.Series(filtered_sWGS_table['Study'])
sample = pd.Series(filtered_sWGS_table['Sample'])
patients = list(set(study_id))

