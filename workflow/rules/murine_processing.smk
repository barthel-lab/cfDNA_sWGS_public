# rule RemoveEmptyReads:
# 	input:
# 		f1 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[5],
# 		f2 = lambda wildcards: filtered_sWGS_table.loc[wildcards.sampleid].iloc[6]
# 	output:
# 		p1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
# 		p2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz",
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/ExtractUmis/{patient_id}/{sampleid}.cutadapt.log"
# 	shell:
# 		"""
# 			cutadapt -m 7 -o {output.p1} -p {output.p2} {input.f1} {input.f2} &> {log}
# 		"""
# # Step 1.2: FASTQ -> uBam 5M (5 bases for UMI) + 2S (2 bases for spacer) T (rest is DNA sequence)
# rule fastq_to_ubam:
# 	input:
# 		r1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
# 		r2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz"
# 	params:
# 		rs1 = "5M2S+T",
# 		rs2 = "5M2S+T",
# 		lane = "L001",
# 		sampleName = "{sampleid}"
# 	output:
# 		bam = base_path + "{study_id}/bam_processing/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam"
# 	resources:
# 		mem_gb = 1
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/fastq_to_ubam/{patient_id}/{sampleid}.log"
# 	shell:
# 		"""
# 			fgbio --compression 1 --async-io FastqToBam \
# 				--input {input.r1} {input.r2} \
# 				--read-structures {params.rs1} {params.rs2} \
# 				--umi-tag RX \
# 				--sample {params.sampleName} \
# 				--library {params.sampleName} \
# 				--platform-unit {params.lane} \
# 				--output {output.bam} &> {log}
# 		"""

# rule align_bam_hg38:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam",
# 		fasta = human
# 	output:
# 		bam = base_path + "{study_id}/xenomapper/align/{patient_id}/{sampleid}_hg38.bam"
# 	threads:
# 		1
# 	resources:
# 		mem_gb = 14
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/align_bam/{patient_id}/{sampleid}_hg38.log"
# 	shell:
# 		"""
# 		samtools fastq {input.bam} \
# 			| bwa mem -t {threads} -p -K 150000000 -Y {input.fasta} - \
# 			| fgbio --compression 1 --async-io ZipperBams \
# 			--unmapped {input.bam} \
# 			--ref {input.fasta} \
# 			--output {output.bam} \
# 			&> {log}
# 		"""

# rule align_bam_mm10:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam",
# 		fasta = mouse
# 	output:
# 		bam = base_path + "{study_id}/xenomapper/align/{patient_id}/{sampleid}_mm10.bam"
# 	threads:
# 		1
# 	resources:
# 		mem_gb = 14
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/align_bam/{patient_id}/{sampleid}_mm10.log"
# 	shell:
# 		"""
# 		samtools fastq {input.bam} \
# 			| bwa mem -t {threads} -p -K 150000000 -Y {input.fasta} - \
# 			| fgbio --compression 1 --async-io ZipperBams \
# 			--unmapped {input.bam} \
# 			--ref {input.fasta} \
# 			--output {output.bam} \
# 			&> {log}
# 		"""

# # rule filter_bam:
# # 	input:
# # 		bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}.{genomes}.bam"
# # 	output:
# # 		bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_{genomes}_filtered.bam"
# # 	threads: 8
# # 	log:
# # 		base_path + "{study_id}/bam_processing/logs/filter_bam/{patient_id}/{sampleid}_{genomes}.log"
# # 	shell:
# # 		"""
# # 		samtools view -@ {threads} -f 2 -F 256 -u {input.bam} | \
# # 		samtools sort -n -@ {threads} -o {output.bam} - &> {log}
# # 		"""

# # rule xenomapper:
# # 	input:
# # 		human = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_hg38_filtered.bam",
# # 		mouse = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_mm10_filtered.bam"
# # 	output:
# # 		human = base_path + "{study_id}/XenoSplit/{patient_id}/{sampleid}_hg38.bam"
# # 	params:
# # 		outdir = base_path + "{study_id}/XenoSplit/{patient_id}/{sampleid}"
# # 	shell:
# # 		"""
# # 		python2 /tgen_labs/barthel/software/XenoSplit/xenosplit.py --pairedEnd --out {params.outdir} {input.human} {input.mouse} 
# # 		"""
# rule xenomapper:
# 	input:
# 		human = "/scratch/smankame/Murine/xenomapper/align/GBMMurine_0009/GBM_VEH_M0_hg38.bam",
# 		mouse = "/scratch/smankame/Murine/xenomapper/align/GBMMurine_0009/GBM_VEH_M0_mm10.bam"
# 	output:
# 		human = "/scratch/smankame/Murine/xenomapper/output_test/GBMMurine_0009/GBM_VEH_M0_primary_specific.bam",
# 		mouse ="/scratch/smankame/Murine/xenomapper/output_test/GBMMurine_0009/GBM_VEH_M0_secondary_specific.bam"
# 	params:
# 		outdir = "/scratch/smankame/Murine/xenomapper/output_test/GBMMurine_0009/GBM_VEH_M0"
# 	shell:
# 		"""
# 		xenomapper2 \
# 			--primary {input.human} \
# 			--secondary {input.mouse} \
# 			--conservative  \
# 			--min-score 30 \
# 			--basename {params.outdir}
# 		"""

# # Step 1.3: uBam -> Mapped BAM
# rule align_bam_combined:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam",
# 		fasta = ref
# 	output:
# 		bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_combined.bam"
# 	threads:
# 		16
# 	resources:
# 		mem_gb = 14
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/align_bam/{patient_id}/{sampleid}_combined.log"
# 	shell:
# 		"""
# 		samtools fastq {input.bam} \
# 			| bwa mem -t {threads} -p -K 150000000 -Y {input.fasta} - \
# 			| fgbio --compression 1 --async-io ZipperBams \
# 			--unmapped {input.bam} \
# 			--ref {input.fasta} \
# 			--output {output.bam} \
# 			&> {log}
# 		"""
# # 

# # rule filter_primary:
# #     input:
# #         bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_combined.bam"
# #     output:
# #         bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_primary.bam"
# #     shell:
# #         """
# #         samtools view -b -F 2304 {input.bam} > {output.bam}
# #         """

# # Step 1.4: Mapped BAM -> Grouped BAM  **deduplication step**
# rule group_reads:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_combined.bam",
# 	output:
# 		bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_grouped.bam",
# 		stats = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_grouped-family-sizes.txt"
# 	params:
# 		allowed_edits = 1,
# 	threads:
# 		2
# 	resources:
# 		mem_gb = 8
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/group_reads/{patient_id}/{sampleid}.log"
# 	shell:
# 		"""
# 			fgbio --compression 1 --async-io GroupReadsByUmi \
# 				--input {input.bam} \
# 				--strategy paired \
# 				--allow-inter-contig false \
# 				--edits {params.allowed_edits} \
# 				--raw-tag RX \
# 				--assign-tag MI \
# 				--min-map-q 20 \
# 				--output {output.bam} \
# 				--family-size-histogram {output.stats} \
# 				&> {log}
# 		"""

# ## Phase 2(a): GroupedBam -> Filtered Consensus

# # Step 2(a).1: GroupedBam -> Consensus uBam
# rule call_consensus_reads:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/align/{patient_id}/{sampleid}_grouped.bam",
# 	output:
# 		bam = base_path + "{study_id}/bam_processing/call_consensus_reads/{patient_id}/{sampleid}.cons.unmapped.bam"
# 	params:
# 		min_reads = 1,
# 		min_base_qual = 10
# 	threads:
# 		4
# 	resources:
# 		mem_gb = 8
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/call_consensus_reads/{patient_id}/{sampleid}.log"
# 	shell:
# 		"""
# 			fgbio --compression 1 CallMolecularConsensusReads \
# 				--input {input.bam} \
# 				--output {output.bam} \
# 				--min-reads {params.min_reads} \
# 				--tag MI \
# 				--error-rate-pre-umi 45 \
# 				--error-rate-post-umi 40 \
# 				--min-input-base-quality {params.min_base_qual} \
# 				--threads {threads} \
# 				&> {log}
# 		"""

# # Step 2(a).2: Consensus uBam -> Consensus Mapped BAM
# rule realign_consensus_reads:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/call_consensus_reads/{patient_id}/{sampleid}.cons.unmapped.bam",
# 		fasta = ref
# 	output:
# 		bam = base_path + "{study_id}/bam_processing/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.bam"
# 	resources:
# 		mem_gb = 4
# 	threads:
# 		16
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.log"
# 	shell:
# 		"""
# 			samtools fastq {input.bam} \
# 				| bwa mem -t 16 -p -K 100000000 -Y {input.fasta} - \
# 				| fgbio --compression 1 --async-io ZipperBams \
# 				--unmapped {input.bam} \
# 				--ref {input.fasta} \
# 				--tags-to-reverse Consensus \
# 				--tags-to-revcomp Consensus \
# 				--output {output.bam} &> {log}
# 		"""

# # Step 2(a).3: Consensus Mapped -> Consensus Filtered & Sorted BAM
# rule sort_consensus_reads:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/realign_consensus_reads/{patient_id}/{sampleid}.cons.mapped.bam",
# 		fasta = ref
# 	output:
# 		bam = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}.bam"
# 	params:
# 		min_reads = 1,
# 		min_base_qual = 40,
# 		max_error_rate = 0.2
# 	threads:
# 		8
# 	resources:
# 		mem_gb = 8
# 	log:
# 		base_path + "{study_id}/bam_processing/logs/filter_consensus_reads/{patient_id}/{sampleid}.log"
# 	shell:
# 		"""
# 			(samtools sort --threads {threads} {input.bam} \
# 				-o {output.bam}) &> {log}
# 		"""

# rule separate_reads:
#     input:
#         bam = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}.bam"
#     output:
#         hg = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}_hg38.bam",
#         mm = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}_mm39.bam",
#         amb = base_path + "{study_id}/bam_processing/sort_consensus_reads/{patient_id}/{sampleid}_subset_amb.bam"
#     threads: 16
#     resources:
#         mem_mb = 150000
#     log:
#         base_path + "{study_id}/bam_processing/logs/sort_consensus_reads/{patient_id}/{sampleid}.log"
#     shell:
#         """
#         python /tgen_labs/barthel/software/github/barthel/cfDNA_sWGS_public/workflow/scripts/separate_reads.py {input.bam} {output.hg} {output.mm} {output.amb} &> {log}
#         """