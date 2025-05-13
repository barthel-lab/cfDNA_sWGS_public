import itertools 
import random 
import pandas as pd

### references needed for pipeline
# hg38 = "/insert/path/to/Homo_sapiens_assembly38.fasta"
# vcf = "/insert/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"
hg38 = "/tgen_labs/barthel/references/GRCh38/hg38_all_viral.fasta"
vcf = "/tgen_labs/barthel/references/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf"
workdir: "/scratch/smankame/slurm/"  ##this is where all the slurm files will be created
hotspot_file = "/tgen_labs/barthel/references/GRCh38/lifted_hotspot_idh1_tert_grc37.vcf"
germline = "/tgen_labs/barthel/Portnow_COH/ref/af-only-gnomad.hg38.vcf.gz"
nPON = "/tgen_labs/barthel/Portnow_COH/ref/xsomatic-hg38_1000g_pon.hg38.vcf"
variant = "/home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz"
intervals = "/tgen_labs/barthel/Portnow_COH/ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
fun_lib = "/tgen_labs/barthel/references/GRCh38/funcotator_dataSources.v1.7.20200521s"


### this path is where all your files will be made
base_path = "/scratch/smankame/"

### user makes a csv file where the first column is study name, the second column is sample name and the third and fourth columns are fastq_R1 and fastq_R2 files.
### study name (first column) will be used to name the folders and sample name (second column) will be used to name all subsequent files
sWGS_table = pd.read_csv("/tgen_labs/barthel/software/github/barthel/cfDNA/sWGS/config/sWGS_data_table.txt",sep='\t', header=None, names=["Study", "Patient", "Sample", "R1","R2", "OldName-R1", "OldName-R2", "Library_Type"])
sWGS_table.index = sWGS_table['Sample']

### if there are multiple studies within the csv table, patient filter will work with specific subsets of samples. Change the study filter value to the patient name. 
study_filter = "Portnow"
filtered_sWGS_table = sWGS_table[sWGS_table['Study'] == study_filter]

study_id=pd.Series(filtered_sWGS_table['Study'])
all_samples = pd.Series(filtered_sWGS_table['Sample'])
patient_id =pd.Series(filtered_sWGS_table['Patient'])


study_list = filtered_sWGS_table['Study'].tolist()
patient_list = filtered_sWGS_table['Patient'].tolist()
sample_list = filtered_sWGS_table['Sample'].tolist()

# print(study_list)
# print(patient_list)
print(sample_list)

# output_files = expand(base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{sampleid}_BQSR_hg38.bam", zip, study_id = study_list, patient_id = patient_list, sampleid=sample_list),
# print("Expanded output files:")
# for f in output_files:
# 	print(f)
# 	print(len(f))