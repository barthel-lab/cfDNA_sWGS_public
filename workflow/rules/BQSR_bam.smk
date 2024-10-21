# Step 1.1: Remove empty reads
rule RemoveEmptyReads:
	input:
		f1 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[2],
		f2 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[3]
	output:
		p1 = base_path + "{study_id}/umi-fgbio/ExtractUmis/{sampleid}/{sampleid}_processed.p1.fastq.gz",
		p2 = base_path + "{study_id}/umi-fgbio/ExtractUmis/{sampleid}/{sampleid}_processed.p2.fastq.gz"
	log:
		base_path + "{study_id}/umi-fgbio/logs/ExtractUmis/{sampleid}.cutadapt.log"
	shell:
		"""
			cutadapt -m 7 -o {output.p1} -p {output.p2} {input.f1} {input.f2} \
				2> {log}
		"""

# Step 1.2: FASTQ -> uBam
rule fastq_to_ubam:
	"""Generates a uBam from R1 and R2 fastq files."""
	input:
		r1 = base_path + "{study_id}/umi-fgbio/ExtractUmis/{sampleid}/{sampleid}_processed.p1.fastq.gz",
		r2 = base_path + "{study_id}/umi-fgbio/ExtractUmis/{sampleid}/{sampleid}_processed.p2.fastq.gz"
	params:
		rs1 = "5M2S+T",
		rs2 = "5M2S+T",
		lane = "L001",
		sampleName = "{sampleid}"
	output:
		bam = base_path + "{study_id}/umi-fgbio/fastq_to_ubam/{sampleid}.unmapped.bam"
	resources:
		mem_gb = 1
	log:
		base_path + "{study_id}/umi-fgbio/logs/fastq_to_ubam/{sampleid}.log"
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
rule align_bam:
	"""Takes an unmapped BAM and generates an aligned BAM using bwa and ZipperBams."""
	input:
		bam = base_path + "{study_id}/umi-fgbio/fastq_to_ubam/{sampleid}.unmapped.bam",
		fasta = hg38
	output:
		bam = base_path + "{study_id}/umi-fgbio/align/{sampleid}.mapped.bam"
	threads:
		16
	resources:
		mem_gb = 14
	log:
		base_path + "{study_id}/umi-fgbio/logs/align_bam/{sampleid}.log"
	shell:
		"""
			(samtools fastq {input.bam} \
				| bwa mem -t {threads} -p -K 150000000 -Y {input.fasta} - \
				| fgbio --compression 1 --async-io ZipperBams \
				--unmapped {input.bam} \
				--ref {input.fasta} \
				--output {output.bam} \
				&> {log})
		"""

# Step 1.4: Mapped BAM -> Grouped BAM
rule group_reads:
	"""Group the raw reads by UMI and position ready for consensus calling."""
	input:
		bam = base_path + "{study_id}/umi-fgbio/align/{sampleid}.mapped.bam",
	output:
		bam = base_path + "{study_id}/umi-fgbio/align/{sampleid}/{sampleid}.grouped.bam",
		stats = base_path + "{study_id}/umi-fgbio/align/{sampleid}/{sampleid}.grouped-family-sizes.txt"
	params:
		allowed_edits = 1,
	threads:
		2
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/umi-fgbio/logs/group_reads/{sampleid}.log"
	shell:
		"""
			fgbio --compression 1 --async-io GroupReadsByUmi \
				--input {input.bam} \
				--strategy paired \
				--allow-inter-contig true \
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
rule call_consensus_reads:
	"""Call consensus reads from the grouped reads."""
	input:
		bam = base_path + "{study_id}/umi-fgbio/align/{sampleid}/{sampleid}.grouped.bam",
	output:
		bam = base_path + "{study_id}/umi-fgbio/call_consensus_reads/{sampleid}/{sampleid}.cons.unmapped.bam"
	params:
		min_reads = 1,
		min_base_qual = 10
	threads:
		4
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/umi-fgbio/logs/call_consensus_reads/{sampleid}.log"
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
rule realign_consensus_reads:
	input:
		bam = base_path + "{study_id}/umi-fgbio/call_consensus_reads/{sampleid}/{sampleid}.cons.unmapped.bam",
		fasta = hg38
	output:
		bam = base_path + "{study_id}/umi-fgbio/realign_consensus_reads/{sampleid}/{sampleid}.cons.mapped.bam"
	resources:
		mem_gb = 4
	threads:
		16
	log:
		base_path + "{study_id}/umi-fgbio/logs/realign_consensus_reads/{sampleid}.log"
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
rule sort_consensus_reads:
	"""Sorts consensus reads into coordinate order."""
	input:
		bam = base_path + "{study_id}/umi-fgbio/realign_consensus_reads/{sampleid}/{sampleid}.cons.mapped.bam",
		fasta = hg38
	output:
		bam = base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam"
	params:
		min_reads = 1,
		min_base_qual = 40,
		max_error_rate = 0.2
	threads:
		8
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/umi-fgbio/logs/filter_consensus_reads/{sampleid}.log"
	shell:
		"""
			(samtools sort --threads {threads} {input.bam} \
				-o {output.bam}) &> {log}
		"""


rule CollectDuplicateMetrics:
	input:
		base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam"
	output:
		base_path + "{study_id}/umi-fgbio/QC/{sampleid}.CollectDuplicateMetrics.mapped._grouped.txt"
	log:
		base_path + "{study_id}/umi-fgbio/logs/QC/mapped_grouped_{sampleid}CollectDuplicateMetrics.log"
	shell:
		"""
			java -jar /tgen_labs/barthel/software/picard_1.8.jar CollectDuplicateMetrics \
				--INPUT {input} \
				--METRICS_FILE {output} > {log} 2>&1
		"""

rule CollectAlignmentSummaryMetrics:
	input:
		bam = base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam",
		fasta = hg38
	output:
		base_path + "{study_id}/umi-fgbio/QC/{sampleid}_CollectAlignmentSummaryMetrics.txt"
	log:
		base_path + "{study_id}/umi-fgbio/logs/QC/{sampleid}_CollectAlignmentSummaryMetrics.log"
	resources:
		mem_mb = 20000
	shell:
		"""
			java -jar /tgen_labs/barthel/software/picard_1.8.jar CollectAlignmentSummaryMetrics \
				-R {input.fasta} \
				-I {input.bam} \
				-O {output} \
				> {log} 2>&1
		"""

rule fastqc:
	input:
		base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam"
	output:
		base_path + "{study_id}/umi-fgbio/QC/{sampleid}/{sampleid}.cons.mapped.sorted_fastqc.html"
	params:
		dir = base_path + "{study_id}/umi-fgbio/QC/{sampleid}"
	log:
		base_path + "{study_id}/umi-fgbio/logs/fastqc/{sampleid}_grouped_fastqc.log"
	shell:
		"""
			#module add fastqc/0.11.8
			fastqc \
				--extract \
				-o {params.dir} \
				-f bam \
				{input} \
				> {log} 2>&1
		"""

rule CollectInsertSizeMetrics:
	input:
		base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam"
	output:
		metrics = base_path + "{study_id}/umi-fgbio/QC/{sampleid}_CollectInsertSizeMetrics.txt",
		pdf = base_path + "{study_id}/umi-fgbio/QC/{sampleid}_CollectInsertSizeMetrics.pdf"
	log:
		base_path + "{study_id}/umi-fgbio/logs/QC/mapped_grouped_{sampleid}_CollectInsertSizeMetrics.log"
	shell:
		"""
			java -jar /tgen_labs/barthel/software/picard_1.8.jar CollectInsertSizeMetrics \
				I={input} O={output.metrics} \
				H={output.pdf} M=0.5 \
				> {log} 2>&1
		"""

rule multiplemetrics:
	input:
		bam = base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam",
		fasta = hg38
	output:
		base_path + "{study_id}/umi-fgbio/QC/{sampleid}.alignment_summary_metrics.txt"
	log:
		base_path + "{study_id}/umi-fgbio/logs/multiplemetrics/{sampleid}.MultipleMetrics.log"
	message:
		"Computing Multiple Metrics\n"
		"Sample: {wildcards.sampleid}"
	shell:
		"""gatk --java-options -Xmx6g CollectMultipleMetrics \
			-R {input.fasta} \
			-I {input.bam} \
			-O {output} \
			> {log} 2>&1"""

rule indexBAM:
	input:
		bam = base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.cons.mapped.sorted.bam"
	output:
		sort = base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.srt.cons.mapped.sorted.bam"
	threads: 8
	shell:
		"""
			samtools sort -@ {threads} {input.bam} -o {output.sort};
			samtools index {output.sort}
		"""

rule CollectWgsMetrics:
	input:
		bam = base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.srt.cons.mapped.sorted.bam",
		fasta = hg38
	output:
		metrics = base_path + "{study_id}/umi-fgbio/QC/{sampleid}_collect_wgs_metrics.txt"
	log:
		base_path + "{study_id}/umi-fgbio/s/DepthOfCoverage/{sampleid}.CollectWgsMetrics.log"
	shell:
		"""
			java -jar /tgen_labs/barthel/software/picard_1.8.jar CollectWgsMetrics \
				-I {input.bam} -O {output.metrics} \
				-R {input.fasta} > {log} 2>&1
		"""

rule addReadGroups:
	input:
		base_path + "{study_id}/umi-fgbio/sort_consensus_reads/{sampleid}/{sampleid}.srt.cons.mapped.sorted.bam"
	output:
		base_path + "{study_id}/umi-fgbio/readGroups/{sampleid}/{sampleid}_RG_hg38.bam"
	log:
		base_path + "{study_id}/umi-fgbio/logs/readGroups/{sampleid}_RG.log"
	message:
		"Adding RG information in sample {wildcards.sampleid}."
	params:
		paramID=r"{sampleid}"
	resources:
		mem_mb = 50000
	shell:
		"""java -jar /tgen_labs/barthel/software/picard_1.8.jar AddOrReplaceReadGroups \
			-I {input} -O {output} \
			-PL ILLUMINA -RGLB {params.paramID} --RGID {params.paramID} \
			-RGPU {params.paramID} -RGSM {params.paramID} \
			> {log} 2>&1"""

# after addReadGroups calculate: CrosscheckFingerprints

rule BQSR:
	input:
		bam = base_path + "{study_id}/umi-fgbio/readGroups/{sampleid}/{sampleid}_RG_hg38.bam",
		fasta = hg38,
		dbsnp = vcf
	output:
		base_path + "{study_id}/umi-fgbio/BQSR/{sampleid}/{sampleid}_recall.table"
	message:
		"Running BaseRecalibrator for {wildcards.sampleid}."
	log:
		base_path + "{study_id}/umi-fgbio/logs/BQSR/{sampleid}_recal.table.log"
	resources:
		mem_mb = 80000
	shell:
		"""gatk BaseRecalibrator \
			--use-original-qualities \
			-I {input.bam} \
			-R {input.fasta} \
			--known-sites {input.dbsnp} \
			-O {output} > {log} 2>&1"""


# gather BQSR reports
rule GatherBQSRReports:
	input:
		base_path + "{study_id}/umi-fgbio/BQSR/{sampleid}/{sampleid}_recall.table"
	output:
		base_path + "{study_id}/umi-fgbio/QC/{sampleid}/{sampleid}_recall.data.table"
	log:
		base_path + "{study_id}/umi-fgbio/logs/GatherBQSRReports/{sampleid}.GatherBQSRReports.log"
	shell:
		"""
			gatk GatherBQSRReports \
				--input {input} \
				--output {output} > {log} 2>&1
		"""

# apply BQSR
rule ApplyBQSR:
	input:
		bam = base_path + "{study_id}/umi-fgbio/readGroups/{sampleid}/{sampleid}_RG_hg38.bam",
		table = base_path + "{study_id}/umi-fgbio/BQSR/{sampleid}/{sampleid}_recall.table",
		fasta = hg38
	output:
		base_path + "{study_id}/umi-fgbio/BQSR/{sampleid}/{sampleid}_BQSR_hg38.bam"
	message:
		"Applying quality filter in {wildcards.sampleid}."
	log:
		base_path + "{study_id}/umi-fgbio/logs/{sampleid}_applyBQSR.bam.log"
	resources:
		mem_mb = 80000
	shell:
		"""gatk ApplyBQSR \
			-R {input.fasta} \
			-I {input.bam} \
			--bqsr-recal-file {input.table} \
			-O {output} > {log} 2>&1"""