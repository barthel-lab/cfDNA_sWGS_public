rule xenome_index:
	input:
		human = human,
		mouse = mouse
	output:
		base_path + "{study_id}/xenome/index"
	log:
		base_path + "{study_id}/xenome/index.log"
	resources:
		runtime = 480,
		mem_mb = 50000
	threads: 8
	shell:
		"/tgen_labs/barthel/software/gossamer/build/src/xenome index -K 31 -T {threads} -P {output} -H {input.mouse} -G {input.human} &> {log}"

rule xenome_RemoveEmptyReads:
	input:
		f1 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[5],
		f2 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[6]
	output:
		p1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		p2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz",
	log:
		base_path + "{study_id}/bam_processing/logs/ExtractUmis/{patient_id}/{sampleid}.cutadapt.log"
	shell:
		"""
			cutadapt -m 7 -o {output.p1} -p {output.p2} {input.f1} {input.f2} &> {log}
		"""

rule xenome_classify:
	input:
		r1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		r2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz"
	params:
		ref =base_path + "{study_id}/xenome/index"
	output:
		base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_"
	log:
		base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_classify.log"
	resources:
		runtime = 480,
		mem_mb = 50000
	threads: 8
	shell:
		"/tgen_labs/barthel/software/gossamer/build/src/xenome classify -P {params.ref} -i {input.r1} -i {input.r2} --pairs --output-filename-prefix {output}"

# Step 1.2: FASTQ -> uBam 5M (5 bases for UMI) + 2S (2 bases for spacer) T (rest is DNA sequence)
rule xenome_fastq_to_ubam:
	input:
		r1 = base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}__host_1.fastq",
		r2 = base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}__host_2.fastq"
	params:
		rs1 = "5M2S+T",
		rs2 = "5M2S+T",
		lane = "L001",
		sampleName = "{sampleid}"
	output:
		bam = base_path + "{study_id}/xenome/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam"
	resources:
		mem_gb = 1
	log:
		base_path + "{study_id}/xenome/logs/fastq_to_ubam/{patient_id}/{sampleid}.log"
	shell:
		"""
			fgbio --compression 1 --async-io FastqToBam \
				--input {input.r1} {input.r2} \
				--read-structures {params.rs1} {params.rs2} \
				--umi-tag RX \
				--sample {params.sampleName} \
				--library {params.sampleName} \
				--platform-unit {params.lane} \
				--output {output.bam} &> {log}
		"""

# Step 1.3: uBam -> Mapped BAM
rule xenome_align_bam:
	input:
		bam = base_path + "{study_id}/xenome/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}_hg38.bam"
	threads:
		16
	resources:
		mem_gb = 14
	log:
		base_path + "{study_id}/xenome/logs/align_bam/{patient_id}/{sampleid}_combined.log"
	shell:
		"""
		samtools fastq {input.bam} \
			| bwa mem -t {threads} -p -K 150000000 -Y {input.fasta} - \
			| fgbio --compression 1 --async-io ZipperBams \
			--unmapped {input.bam} \
			--ref {input.fasta} \
			--output {output.bam} \
			&> {log}
		"""

# Step 1.4: Mapped BAM -> Grouped BAM  **deduplication step**
rule xenome_group_reads:
	input:
		bam = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}_hg38.bam",
	output:
		bam = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}_grouped.bam",
		stats = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}_grouped-family-sizes.txt"
	params:
		allowed_edits = 1,
	threads:
		2
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/xenome/logs/group_reads/{patient_id}/{sampleid}.log"
	shell:
		"""
			fgbio --compression 1 --async-io GroupReadsByUmi \
				--input {input.bam} \
				--strategy paired \
				--allow-inter-contig false \
				--edits {params.allowed_edits} \
				--raw-tag RX \
				--assign-tag MI \
				--min-map-q 20 \
				--output {output.bam} \
				--family-size-histogram {output.stats} \
				&> {log}
		"""

## Phase 2(a): GroupedBam -> Filtered Consensus

# Step 2(a).1: GroupedBam -> Consensus uBam
rule xenome_all_consensus_reads:
	input:
		bam = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}_grouped.bam",
	output:
		bam = base_path + "{study_id}/xenome/call_consensus_reads/{patient_id}/{sampleid}.cons.unmapped.bam"
	params:
		min_reads = 1,
		min_base_qual = 10
	threads:
		4
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/xenome/logs/call_consensus_reads/{patient_id}/{sampleid}.log"
	shell:
		"""
			fgbio --compression 1 CallMolecularConsensusReads \
				--input {input.bam} \
				--output {output.bam} \
				--min-reads {params.min_reads} \
				--tag MI \
				--error-rate-pre-umi 45 \
				--error-rate-post-umi 40 \
				--min-input-base-quality {params.min_base_qual} \
				--threads {threads} \
				&> {log}
		"""

# Step 2(a).2: Consensus uBam -> Consensus Mapped BAM
rule xenome_realign_consensus_reads:
	input:
		bam = base_path + "{study_id}/xenome/call_consensus_reads/{patient_id}/{sampleid}.cons.unmapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/xenome/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.bam"
	resources:
		mem_gb = 4
	threads:
		16
	log:
		base_path + "{study_id}/xenome/logs/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.log"
	shell:
		"""
			samtools fastq {input.bam} \
				| bwa mem -t 16 -p -K 100000000 -Y {input.fasta} - \
				| fgbio --compression 1 --async-io ZipperBams \
				--unmapped {input.bam} \
				--ref {input.fasta} \
				--tags-to-reverse Consensus \
				--tags-to-revcomp Consensus \
				--output {output.bam} &> {log}
		"""

# Step 2(a).3: Consensus Mapped -> Consensus Filtered & Sorted BAM
rule xenome_sort_consensus_reads:
	input:
		bam = base_path + "{study_id}/xenome/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/xenome/sort_consensus_reads/{patient_id}/{sampleid}.bam"
	params:
		min_reads = 1,
		min_base_qual = 40,
		max_error_rate = 0.2
	threads:
		8
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/xenome/logs/filter_consensus_reads/{patient_id}/{sampleid}.log"
	shell:
		"""
			(samtools sort --threads {threads} {input.bam} \
				-o {output.bam}) &> {log}
		"""

###run CNA analysis to see if there's any GBM related alternation and also run mutect2 to see which samples show the GBM43 variants

