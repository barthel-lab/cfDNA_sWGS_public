# cfDNA public repo

From sWGS sequencing data, this pipeline includes:

The hg38 fasta file (Homo_sapiens_assembly38.fasta) needed for alignment is available through this link: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?inv=1&invt=AbxUZg&prefix=&forceOnObjectsSortingFiltering=false

1.  Creating BQSR bam files and generating QC reports from raw sample sequencing fastq files
Preprocessing:
  -rename_fastq: Rename and organize input FASTQs
  - RemoveEmptyReads: Remove reads < 7bp using cutadapt
UMI Extraction & Mapping:
  - fastq_to_ubam: Convert to uBam with read structure 5M2S+T (UMI+spacer+read)
  - align_bam: Map reads to hg38 and stitch in with ZipperBams
  - group_reads: UMI-aware grouping for consensus calling
Consensus Calling:
  - call_consensus_reads: Call consensus with fgbio
  - realign_consensus_reads: Re-align consensus reads
  - sort_consensus_reads: Sort mapped consensus reads
Read Group Assignment:
  - addReadGroups: Add RG (read group) info to mapped consensus BAMs using Picard.
Base Quality Score Recalibration (BQSR):
  - BQSR: Generate recalibration table using GATK BaseRecalibrator with known variants (dbSNP).
  - ApplyBQSR: Apply recalibration to BAM file using the generated recal table.
  - BQSR_index: Index the recalibrated BAM using samtools.

2.  Generating ichorCNA files for copy number analysis and tumor fraction quantification
  - If patient normal is not available, use PON file from ichorCNA
  - directions to download ichorCNA are provided at: https://github.com/broadinstitute/ichorCNA/wiki/Installation
  - file paths to downloaded ichorCNA bin folder will need to be updated within ichorCNA.smk

3. Running variant calling on WGS data in multisample mode
  - reference files for this pipeline are available through the following links:
    
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?inv=1&invt=AbxUZg&prefix=&forceOnObjectsSortingFiltering=false
 -  resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
 - Homo_sapiens_assembly38.dbsnp138.vcf
 
 https://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/
 - gnomad.genomes.r3.0.sites.vcf.gz     
 
 https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?inv=1&invt=AbxUcQ&prefix=&forceOnObjectsSortingFiltering=false
 - af-only-gnomad.hg38.vcf.gz
 - 1000g_pon.hg38.vcf.gz