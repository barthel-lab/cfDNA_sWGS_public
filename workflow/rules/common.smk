import itertools 
import random 
import pandas as pd

username = "wsmith" ###change based on user
workdir: "/scratch/{username}/slurm/"  ##this is where all the slurm files will be created

### this path is where all your files will be made
base_path = "/scratch/{username}/"
picard_path = "/tgen_labs/barthel/software/picard_1.8.jar"
labs_dir ="/tgen_labs/barthel/projects/"

####references for BQSR pipeline
hg38 = "/tgen_labs/barthel/references/GRCh38/hg38_all_viral.fasta"
ref = "/tgen_labs/barthel/projects/GBM_Cell_Culture/ref/hg38_mm39.fasta"

####references for variant calling pipeline
vcf = "/tgen_labs/barthel/references/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf"
hotspot_file = "/tgen_labs/barthel/references/GRCh38/lifted_hotspot_idh1_tert_grc37.vcf"
germline = "/tgen_labs/barthel/projects/Portnow_COH/ref/af-only-gnomad.hg38.vcf.gz"
nPON = "/tgen_labs/barthel/projects/Portnow_COH/ref/xsomatic-hg38_1000g_pon.hg38.vcf"
variant = "/home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz"
intervals = "/tgen_labs/barthel/projects/Portnow_COH/ref/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
fun_lib = "/tgen_labs/barthel/references/GRCh38/funcotator_dataSources.v1.7.20200521s"

####references/filepaths for ichorCNA pipeline
ichorcna_path = "/home/{username}/miniforge3/envs/ichorcna/bin/ichorCNA" 
readCounter_path = "/home/{username}/miniforge3/envs/ichorcna/bin/readCounter"
bam_file_path = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}.srt" ##connect the bam files created in BQSR pipeline straight to ichorCNA pipeline
vcf_intervals = "/tgen_labs/barthel/references/GRCh38/1000G_phase1.snps.high_confidence.hg38.vcf.interval_list"
gatk_intervals_1000bp = "/tgen_labs/barthel/references/GRCh38/Homo_sapiens_assembly38.1000bp.interval_list"
gatk_intervals_100kbp = "/tgen_labs/barthel/references/GRCh38/Homo_sapiens_assembly38.100kbp.interval_list" ##custom interval list can be made using GATK PreprocessIntervals (rule preprocessIntervals in ichorCNA.smk needs human genome fasta file)
hg38_dict = "/tgen_labs/barthel/references/GRCh38/Homo_sapiens_assembly38.dict"

####references for murine pipeline
mouse = "" ##insert file path for mouse genome fasta file
human = "" ##insert file path for human genome fasta file


### user makes a csv file where the first column is study name, the second column is sample name and the third and fourth columns are fastq_R1 and fastq_R2 files.
### study name (first column) will be used to name the folders and sample name (second column) will be used to name all subsequent files
sWGS_table = pd.read_csv("/tgen_labs/barthel/projects/cfDNA_sWGS_public/config/sWGS_data_table.txt",sep='\t', header=None, names=["Study", "Patient", "Sample", "R1","R2"])
sWGS_table.index = sWGS_table['Sample']

### if there are multiple studies within the csv table, patient filter will work with specific subsets of samples. Change the study filter value to the patient name. 
study_filter = "Murine"
filtered_sWGS_table = sWGS_table[sWGS_table['Study'] == study_filter]

study_id=pd.Series(filtered_sWGS_table['Study'])
all_samples = pd.Series(filtered_sWGS_table['Sample'])
patient_id =pd.Series(filtered_sWGS_table['Patient'])


study_list = filtered_sWGS_table['Study'].tolist()
patient_list = filtered_sWGS_table['Patient'].tolist()
sample_list = filtered_sWGS_table['Sample'].tolist()

genomes = ["human","mouse"]
