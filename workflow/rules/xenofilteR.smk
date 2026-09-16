rule RemoveEmptyReads:
	input:
		f1 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[5],
		f2 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[6]
	output:
		p1 = base_path + "{study_id}/xenofilteR/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		p2 = base_path + "{study_id}/xenofilteR/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz",
	log:
		base_path + "{study_id}/xenofilteR/logs/ExtractUmis/{patient_id}/{sampleid}.cutadapt.log"
	shell:
		"""
			cutadapt -m 7 -o {output.p1} -p {output.p2} {input.f1} {input.f2} &> {log}
		"""

rule align_bam_hg38:
	input:
		p1 = base_path + "{study_id}/xenofilteR/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		p2 = base_path + "{study_id}/xenofilteR/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz",
		fasta = human
	output:
		bam = base_path + "{study_id}/xenofilteR/align/{patient_id}/{sampleid}_hg38.bam"
	threads: 16
	resources:
		mem_gb = 14
	log:
		base_path + "{study_id}/xenofilteR/logs/align_bam/{patient_id}/{sampleid}_hg38.log"
	shell:
		"""
		bwa mem -t {threads} -Y -M {input.fasta} {input.p1} {input.p2} \
		| samtools view -b - \
		| samtools sort -n -@ {threads} -o {output.bam}
		&> {log}
		"""
rule align_bam_mm10:
	input:
		p1 = base_path + "{study_id}/xenofilteR/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		p2 = base_path + "{study_id}/xenofilteR/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz",
		fasta = mouse
	output:
		bam = base_path + "{study_id}/xenofilteR/align/{patient_id}/{sampleid}_mm10.bam"
	threads: 16
	resources:
		mem_gb = 14
	log:
		base_path + "{study_id}/xenofilteR/logs/align_bam/{patient_id}/{sampleid}_mm10.log"
	shell:
		"""
		bwa mem -t {threads} -Y -M {input.fasta} {input.p1} {input.p2} \
		| samtools view -b - \
		| samtools sort -n -@ {threads} -o {output.bam}
		&> {log}
		"""

rule filter_bam:
	input:
		bam = base_path + "{study_id}/xenofilteR/align/{patient_id}/{sampleid}_{genomes}.bam"
	output:
		bam = base_path + "{study_id}/xenofilteR/align/{patient_id}/{sampleid}_{genomes}_filtered.bam"
	threads: 8
	log:
		base_path + "{study_id}/xenofilteR/logs/filter_bam/{patient_id}/{sampleid}_{genomes}.log"
	shell:
		"""
		samtools view -u -b -F 2304 {input.bam} | \
		samtools sort -n -@ {threads} -o {output.bam} - &> {log}
		"""

rule xenofilter:
	input:
		human = base_path + "{study_id}/xenofilteR/align/{patient_id}/{sampleid}_human_filtered.bam",
		mouse = base_path + "{study_id}/xenofilteR/align/{patient_id}/{sampleid}_mouse_filtered.bam"
	output:
        filtered_bam = base_path + "{study_id}/xenofilteR/filtered/{patient_id}/{sampleid}_Filtered.bam"
    params:
        outdir = base_path + "{study_id}/xenofilteR/filtered/{patient_id}/"
    threads: 4
    resources:
        mem_gb = 16
    log:
        base_path + "{study_id}/xenofilteR/logs/xenofilter/{patient_id}/{sampleid}.log"
    script:
        "../scripts/xenofilter.R"

rule bam_to_fastq:
	input:
		base_path + "{study_id}/xenofilteR/filtered/{patient_id}/{sampleid}_Filtered.bam"
	output:
		r1 = base_path + "{study_id}/disambiguate/{patient_id}/{sampleid}_R1.fastq",
		r2 = base_path + "{study_id}/disambiguate/{patient_id}/{sampleid}_R2.fastq"
	params:
		sorted_bam = base_path + "{study_id}/disambiguate/{patient_id}/{sampleid}_sorted.bam"
	resources:
		mem_mb = 15000
	shell:
		"bedtools bamtofastq -i {params.sorted_bam} -fq {output.r1} -fq2 {output.r2}"

rule merge_human_fastqs:
	input:
		xf_R1 = base_path + "{study_id}/xenofilteR/{patient_id}/{sampleid}_R1.fastq",
		xf_R2 = base_path + "{study_id}/xenofilteR/{patient_id}/{sampleid}_R2.fastq",
		xenome_R1 = base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_graft_1.fastq",
		xenome_R2 = base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_graft_2.fastq"
	output:
		R1 = base_path + "{study_id}/merged_human/{patient_id}/{sampleid}_R1.fastq.gz",
		R2 = base_path + "{study_id}/merged_human/{patient_id}/{sampleid}_R2.fastq.gz"
	threads: 2
	shell:
		"""
		cat {input.xf_R1} {input.xenome_R1} \
		  | seqkit rmdup -n \
		  | gzip > {output.R1}

		cat {input.xf_R2} {input.xenome_R2} \
		  | seqkit rmdup -n \
		  | gzip > {output.R2}
		"""

rule fastq_to_ubam:
	input:
		r1 = base_path + "{study_id}/merged_human/{patient_id}/{sampleid}_R1.fastq.gz",
		r2 = base_path + "{study_id}/merged_human/{patient_id}/{sampleid}_R2.fastq.gz"
	params:
		rs1 = "5M2S+T",
		rs2 = "5M2S+T",
		lane = "L001",
		sampleName = "{sampleid}"
	output:
		bam = base_path + "{study_id}/merged_human/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam"
	resources:
		mem_gb = 1
	log:
		base_path + "{study_id}/merged_human/logs/fastq_to_ubam/{patient_id}/{sampleid}.log"
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
# # Step 1.3: uBam -> Mapped BAM
rule align:
	input:
		bam = base_path + "{study_id}/merged_human/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/merged_human/align/{patient_id}/{sampleid}.bam"
	threads:
		16
	resources:
		mem_gb = 14
	log:
		base_path + "{study_id}/merged_human/logs/align_bam/{patient_id}/{sampleid}.log"
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
rule group_reads:
	input:
		bam = base_path + "{study_id}/merged_human/align/{patient_id}/{sampleid}.bam",
	output:
		bam = base_path + "{study_id}/merged_human/align/{patient_id}/{sampleid}_grouped.bam",
		stats = base_path + "{study_id}/merged_human/align/{patient_id}/{sampleid}_grouped-family-sizes.txt"
	params:
		allowed_edits = 1,
	threads:
		2
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/merged_human/logs/group_reads/{patient_id}/{sampleid}.log"
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
rule call_consensus_reads:
	input:
		bam = base_path + "{study_id}/merged_human/align/{patient_id}/{sampleid}_grouped.bam",
	output:
		bam = base_path + "{study_id}/merged_human/call_consensus_reads/{patient_id}/{sampleid}.cons.unmapped.bam"
	params:
		min_reads = 1,
		min_base_qual = 10
	threads:
		4
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/merged_human/logs/call_consensus_reads/{patient_id}/{sampleid}.log"
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
		bam = base_path + "{study_id}/merged_human/call_consensus_reads/{patient_id}/{sampleid}.cons.unmapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/merged_human/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.bam"
	resources:
		mem_gb = 4
	threads:
		16
	log:
		base_path + "{study_id}/merged_human/logs/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.log"
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
	input:
		bam = base_path + "{study_id}/merged_human/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/merged_human/sort_consensus_reads/{patient_id}/{sampleid}.bam"
	params:
		min_reads = 1,
		min_base_qual = 40,
		max_error_rate = 0.2
	threads:
		8
	resources:
		mem_gb = 8
	log:
		base_path + "{study_id}/merged_human/logs/filter_consensus_reads/{patient_id}/{sampleid}.log"
	shell:
		"""
			(samtools sort --threads {threads} {input.bam} \
				-o {output.bam}) &> {log}
		"""