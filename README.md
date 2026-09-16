# cfDNA sWGS Processing Pipeline

This Snakemake pipeline processes shallow whole-genome sequencing (sWGS) data from cell-free DNA (cfDNA). It supports four modular sub-pipelines that can be run independently or together depending on your sample type.

## Pipelines overview

**1. BAM processing** (`rules/BQSR_bam.smk`)
Generates consensus BAM files from raw FASTQ input using UMI-aware deduplication, followed by QC.
- Trim short reads (cutadapt)
- FASTQ → unmapped BAM with UMI extraction (fgbio)
- Align to hg38 (bwa mem + fgbio ZipperBams)
- UMI-based grouping and consensus calling (fgbio)
- Re-alignment and sorting of consensus reads (samtools)
- QC: CollectDuplicateMetrics, CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, FastQC, CollectWgsMetrics, CollectQualityYieldMetrics, MultiQC
- *Optional (for variant calling):* AddReadGroups, BaseRecalibration (BQSR), ApplyBQSR — uncomment these rules if needed

**2. ichorCNA copy number analysis** (`rules/ichorCNA.smk`)
Estimates copy number alterations and tumor fraction from sWGS BAMs.
- Generate read depth wiggle files (readCounter)
- Build a panel of normals from matched PBMCs (optional)
- Run ichorCNA
- Arm-level aneuploidy summary per patient (`scripts/aneuploidy_ichorcna.R`)
- *Note: ichorCNA has not been optimized for cfDNA by default — see [ichorCNA docs](https://github.com/broadinstitute/ichorCNA/wiki) for cfDNA-specific parameter guidance*

**3. Somatic variant calling** (`rules/variant_calling.smk`)
Calls somatic mutations using GATK Mutect2 in single-sample or multi-sample mode.
- Mutect2 → LearnReadOrientationModel → GetPileupSummaries → CalculateContamination → FilterMutectCalls → Funcotator annotation → VAF extraction
- *Note: Mutect2 is optimized for WGS tumor tissue; results on cfDNA may show higher noise. Multi-sample mode is implemented but currently commented out.*

**4. Murine / PDX de-mixing** (`rules/xenome.smk`, `rules/xenofilteR.smk`, `rules/murine_processing.smk`)
Separates human and mouse reads from patient-derived xenograft (PDX) samples. Three approaches are included:
- **Xenome** (k-mer based classification using gossamer)
- **xengsort** (k-mer based, faster alternative)
- **BBSplit** (alignment-based)
- After classification, human reads proceed through the same UMI consensus pipeline as above
- `scripts/separate_reads.py` — BAM-level read separator using chromosome name patterns and XA tags

---

## Requirements

### Software
Install the following tools and ensure they are available on your `PATH`:

| Tool | Purpose |
|------|---------|
| [Snakemake](https://snakemake.readthedocs.io) ≥ 7.0 | Workflow manager |
| [fgbio](https://github.com/fulcrumgenomics/fgbio) | UMI extraction and consensus calling |
| [bwa](https://github.com/lh3/bwa) | Read alignment |
| [samtools](https://www.htslib.org/) | BAM manipulation |
| [cutadapt](https://cutadapt.readthedocs.io) | Adapter/short read trimming |
| [GATK](https://gatk.broadinstitute.org/) ≥ 4.0 | Variant calling, BQSR, QC metrics |
| [Picard](https://broadinstitute.github.io/picard/) | QC metrics (CollectDuplicateMetrics, etc.) |
| [bcftools](https://samtools.github.io/bcftools/) | VCF filtering |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-sample QC |
| [MultiQC](https://multiqc.info/) | Aggregated QC report |
| [ichorCNA](https://github.com/broadinstitute/ichorCNA) | Copy number and tumor fraction |
| R ≥ 4.0 | Aneuploidy analysis script |
| Python ≥ 3.8 with `pandas`, `pysam` | Workflow configuration and scripts |

For the murine pipeline only:
| [xenome / gossamer](https://github.com/data61/gossamer) | k-mer classification |
| [xengsort](https://gitlab.com/genomeinformatics/xengsort) | k-mer classification (alternative) |
| [BBTools / BBSplit](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | Alignment-based classification |
| [KMC](https://github.com/refresh-bio/KMC) | k-mer counting |
| [seqtk](https://github.com/lh3/seqtk) | FASTQ to FASTA conversion |

---

## Reference files

### Included in this repository (`workflow/references/`)
- `Homo_sapiens_assembly38.1000bp.interval_list` — custom GATK interval list (1 kb bins)
- `Homo_sapiens_assembly38.100kbp.interval_list` — custom GATK interval list (100 kb bins)
- `Homo_sapiens_assembly38.dict` — hg38 sequence dictionary
- `lifted_hotspot_idh1_tert_grc37.vcf` — IDH1/TERT hotspot VCF (lifted to hg38)
- `gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz` — gnomAD subset for contamination estimation

### Download separately
Update the paths in `config/config.yaml` (see **Configuration** below) after downloading:

**Broad Institute hg38 bundle** — [Google Cloud bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0):
- `Homo_sapiens_assembly38.fasta` (and `.fai`, `.64.amb`, etc.) — main hg38 reference
- `Homo_sapiens_assembly38.dbsnp138.vcf` — known variants for BQSR
- `resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list` — calling regions

**GATK best practices somatic hg38** — [Google Cloud bucket](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38):
- `af-only-gnomad.hg38.vcf.gz` — germline resource for Mutect2
- `1000g_pon.hg38.vcf.gz` — panel of normals for Mutect2

**gnomAD** — [UCSC Euro mirror](https://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/):
- `gnomad.genomes.r3.0.sites.vcf.gz` (full sites VCF, if needed)

**Funcotator data sources** — download via GATK:
```bash
gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
```

**ichorCNA reference files** — included with the ichorCNA R package (`inst/extdata/`):
- `gc_hg38_1000kb.wig`
- `map_hg38_1000kb.wig`
- `GRCh38.GCA_000001405.2_centromere_acen.txt`
- Built-in PON: `HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds` (use this if no matched PBMC normal is available)

**For PDX / murine samples:** you will also need hg38 and mm39 (GRCm39) FASTA files to build the xenome/xengsort index.

---

## Configuration

### 1. Edit `workflow/rules/common.smk`
All file paths and user-specific settings are set at the top of this file. Update the following variables before running:

```python
# Your HPC username
username = "your_username"

# Where all output files will be written (your scratch or project directory)
base_path = "/scratch/your_username/"

# Path to Picard jar
picard_path = "/path/to/picard.jar"

# hg38 reference FASTA (see Reference files section for download)
hg38 = "/path/to/Homo_sapiens_assembly38.fasta"

# For variant calling
vcf       = "/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"
germline  = "/path/to/af-only-gnomad.hg38.vcf.gz"
nPON      = "/path/to/1000g_pon.hg38.vcf.gz"
variant   = "/path/to/gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz"
intervals = "/path/to/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
fun_lib   = "/path/to/funcotator_dataSources/"

# For ichorCNA
ichorcna_path    = "/path/to/ichorCNA"
readCounter_path = "/path/to/readCounter"

# For murine/PDX pipeline (leave empty if not used)
mouse = ""   # path to mouse genome FASTA
human = ""   # path to human genome FASTA (can be same as hg38 above)
```

### 2. Fill in `config/sWGS_data_table.txt`
This is the sample manifest. Add one row per sample with tab-separated columns — the header row is already included:

| Column | Description |
|--------|-------------|
| `Study` | Study or cohort name — used for top-level output folder naming |
| `Patient` | Patient identifier — used for per-patient sub-folders |
| `Sample` | Sample identifier — used to name all output files |
| `R1` | Absolute path to R1 FASTQ (gzipped) |
| `R2` | Absolute path to R2 FASTQ (gzipped) |


To run on a subset of samples, set `study_filter` in `common.smk` to the `Study` value you
---

## Running the pipeline

### Select which sub-pipeline to run
In `workflow/Snakefile`, uncomment the output targets in `rule all` that correspond to the pipeline you want to run. Each section is labeled (BAM processing, ichorCNA, variant calling, murine). Only one set of targets needs to be active at a time.

### Dry run (recommended first)
```bash
snakemake -n --snakefile workflow/Snakefile
```

### Local run
```bash
snakemake --cores 16 --snakefile workflow/Snakefile
```

### SLURM cluster
```bash
snakemake --snakefile workflow/Snakefile \
  --executor slurm \
  --default-resources slurm_account=<your_account> runtime=480 mem_mb=8000 \
  --jobs 50
```

---

## Output structure

All outputs are written under `base_path/{study_id}/`: