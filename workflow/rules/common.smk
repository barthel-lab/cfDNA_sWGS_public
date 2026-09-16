import itertools 
import random 
import pandas as pd

username = "your_username"  # replace with your HPC username 
workdir: "/scratch/{username}/slurm/"  ##this is where all the slurm files will be created

### this path is where all your files will be made
workdir: "/scratch/{username}/slurm/"   # update to your HPC scratch path
base_path = "/scratch/{username}/"       # update to where you want outputs
picard_path = "/path/to/picard.jar"      # path to your picard jar


####references for BQSR pipeline
hg38 = "/path/to/Homo_sapiens_assembly38.fasta"  # standard hg38 FASTA

####references for murine pipeline
mouse = "" ##insert file path for mouse genome fasta file
human = hg38

####references for variant calling pipeline
vcf          = "/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"
hotspot_file = "workflow/references/lifted_hotspot_idh1_tert_grc37.vcf"  # included in repo
germline     = "/path/to/af-only-gnomad.hg38.vcf.gz"
nPON         = "/path/to/1000g_pon.hg38.vcf.gz"
variant      = "workflow/references/gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz"  # included in repo
intervals    = "/path/to/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
fun_lib      = "/path/to/funcotator_dataSources/"

####references/filepaths for ichorCNA pipeline
ichorcna_path    = "/path/to/ichorCNA"         # wherever ichorCNA is installed
readCounter_path = "/path/to/readCounter"
bam_file_path = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}.srt" ##connect the bam files created in BQSR pipeline straight to ichorCNA pipeline
gatk_intervals_1000bp = "workflow/references/Homo_sapiens_assembly38.1000bp.interval_list"  # included in repo
gatk_intervals_100kbp = "workflow/references/Homo_sapiens_assembly38.100kbp.interval_list"  # included in repo
hg38_dict             = "workflow/references/Homo_sapiens_assembly38.dict"                  # included in repo
vcf_intervals         = "/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.interval_list"


### user makes a csv file where the first column is study name, the second column is sample name and the third and fourth columns are fastq_R1 and fastq_R2 files.
### study name (first column) will be used to name the folders and sample name (second column) will be used to name all subsequent files
sWGS_table = pd.read_csv("config/sWGS_data_table.txt",sep='\t', header=None, names=["Study", "Patient", "Sample", "R1","R2"])
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
