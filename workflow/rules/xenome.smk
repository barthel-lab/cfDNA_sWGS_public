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

# rule xenome_RemoveEmptyReads:
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

# rule xenome_classify:
# 	input:
# 		r1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
# 		r2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz"
# 	params:
# 		ref =base_path + "{study_id}/xenome/index/index"
# 	output:
# 		base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}"
# 	log:
# 		base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_classify.log"
# 	resources:
# 		runtime = 480,
# 		mem_mb = 50000
# 	threads: 8
# 	shell:
# 		"/tgen_labs/barthel/software/gossamer/build/src/xenome classify -P {params.ref} -i {input.r1} -i {input.r2} --pairs --output-filename-prefix {output}"

rule xfastq_to_ubam:
	input:
		r1 = base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_graft_1.fastq.gz",
		r2 = base_path + "{study_id}/xenome/classify/{patient_id}/{sampleid}_graft_2.fastq.gz"
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
# # Step 1.3: uBam -> Mapped BAM
rule xalign:
	input:
		bam = base_path + "{study_id}/xenome/fastq_to_ubam/{patient_id}/{sampleid}.unmapped.bam",
		fasta = human
	output:
		bam = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}.bam"
	threads:
		16
	resources:
		mem_gb = 14
	log:
		base_path + "{study_id}/xenome/logs/align_bam/{patient_id}/{sampleid}.log"
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
rule xgroup_reads:
	input:
		bam = base_path + "{study_id}/xenome/align/{patient_id}/{sampleid}.bam",
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
rule xcall_consensus_reads:
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
rule xrealign_consensus_reads:
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
rule xsort_consensus_reads:
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
rule xengsort:
	input:
		human = human,
		mouse = mouse
	output:
		indexdir = directory("/scratch/smankame/{study_id}/xengsort/")
	log:
		base_path + "{study_id}/xengsort/index.log"
	resources:
		runtime = 720,
		mem_mb = 120000
	shell:
		"xengsort index --index {output.indexdir}index -H {input.mouse} -G {input.human} -k 25 -n 4_500_000_000  &> {log}"

rule xengsort_classify:
	input:
		r1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		r2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz"
	params:
		ref =base_path + "{study_id}/xengsort/xengsortindex",
		outdir = base_path + "{study_id}/xengsort/classify/{patient_id}/{sampleid}"
	output:
		base_path + "{study_id}/xengsort/classify/{patient_id}/{sampleid}.host.1.fq.gz"
	log:
		base_path + "{study_id}/xenogsort/classify/{patient_id}/{sampleid}_classify.log"
	resources:
		runtime = 480,
		mem_mb = 50000
	threads: 8
	shell:
		"xengsort classify --index {params.ref} --fastq {input.r1} --pairs {input.r2} --prefix {params.outdir} --mode coverage &> {log}"

rule bbsplit:
	input:
		human = human,
		mouse = mouse,
		r1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		r2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz"
	output:
		hg38_1 = base_path + "{study_id}/bbsplit/{patient_id}/{sampleid}_hg38_1.fq",
		hg38_2 = base_path + "{study_id}/bbsplit/{patient_id}/{sampleid}_hg38_2.fq",
		mm39_1 = base_path + "{study_id}/bbsplit/{patient_id}/{sampleid}_mm39_1.fq",
		mm39_2 = base_path + "{study_id}/bbsplit/{patient_id}/{sampleid}_mm39_2.fq"
	log:
		base_path + "{study_id}/bbsplit/{patient_id}/{sampleid}.log"
	resources:
		runtime = 480,
		mem_mb = 50000
	threads: 8
	shell:
		"""
		mkdir -p {base_path}{wildcards.study_id}/bbsplit/{wildcards.patient_id}
		/home/smankame/miniforge3/envs/xenofilter2/bin/bbsplit.sh \
			ref={input.human},{input.mouse} \
			in={input.r1} in2={input.r2} \
			basename={base_path}{wildcards.study_id}/bbsplit/{wildcards.patient_id}/{wildcards.sampleid}_%_#.fq \
			ambig2=best	&> {log}
		"""

###### testing out KMC
rule count_kmers_genomes:
	input:
		human_genome = human,
		mouse_genome = mouse,
	output:
		human_kmc =  base_path + "{study_id}/count_kmers/human_all_kmers.kmc_suf",
		mouse_kmc =  base_path + "{study_id}/count_kmers/mouse_all_kmers.kmc_suf"
	params:
		human_kmc =  base_path + "{study_id}/count_kmers/human_all_kmers",
		mouse_kmc = base_path + "{study_id}/count_kmers/mouse_all_kmers",
		temp = base_path + "{study_id}/count_kmers"
	log:
		 base_path + "{study_id}/count_kmers/count_kmers.log"
	resources:
		runtime = 480,
		mem_mb = 50000
	threads: 2
	shell:
		"""
		mkdir -p {params.temp}

		kmc -k31 -ci1 -t{threads} -fm \
			{input.human_genome} \
			{params.human_kmc} \
			{params.temp} \
		&& \
		kmc -k31 -ci1 -t{threads} -fm \
			{input.mouse_genome} \
			{params.mouse_kmc} \
			{params.temp}
		&>> {log}
		"""

rule unique_kmers_both:
	input:
		human_kmc = base_path + "{study_id}/count_kmers/human_all_kmers.kmc_suf",
		mouse_kmc = base_path + "{study_id}/count_kmers/mouse_all_kmers.kmc_suf"
	output:
		human_kmc = base_path + "{study_id}/count_kmers/human_unique_kmers.kmc_suf",
		mouse_kmc = base_path + "{study_id}/count_kmers/mouse_unique_kmers.kmc_suf",
		human_txt = base_path + "{study_id}/count_kmers/human_unique_kmers.txt",
		mouse_txt = base_path + "{study_id}/count_kmers/mouse_unique_kmers.txt"
	params:
		human = base_path + "{study_id}/count_kmers/human_all_kmers",
		mouse = base_path + "{study_id}/count_kmers/mouse_all_kmers",
		human_out = base_path + "{study_id}/count_kmers/human_unique_kmers",
		mouse_out = base_path + "{study_id}/count_kmers/mouse_unique_kmers"
	log:
		base_path + "{study_id}/count_kmers/unique_kmers.log"
	threads: 16
	resources:
		runtime = 480,
		mem_mb = 50000
	shell:
		"""
		kmc_tools simple {params.human} {params.mouse} kmers_subtract {params.human_out} \
		&& kmc_tools transform {params.human_out} dump {output.human_txt} \
		\
		&& kmc_tools simple {params.mouse} {params.human} kmers_subtract {params.mouse_out} \
		&& kmc_tools transform {params.mouse_out} dump {output.mouse_txt}
		"""

rule shared_kmers_human_mouse:
	input:
		human_kmc = base_path + "{study_id}/count_kmers/human_all_kmers.kmc_suf",
		mouse_kmc = base_path + "{study_id}/count_kmers/mouse_all_kmers.kmc_suf"
	output:
		kmc  = base_path + "{study_id}/count_kmers/human_mouse_shared_kmers.kmc_suf",
		text = base_path + "{study_id}/count_kmers/human_mouse_shared_kmers.txt"
	params:
		human = base_path + "{study_id}/count_kmers/human_all_kmers",
		mouse = base_path + "{study_id}/count_kmers/mouse_all_kmers",
		out = base_path + "{study_id}/count_kmers/human_mouse_shared_kmers"
	log:
		base_path + "{study_id}/count_kmers/shared_kmers.log"
	threads: 8
	resources:
		runtime = 240,
		mem_mb = 30000
	shell:
		"""
		kmc_tools simple {params.human} {params.mouse} intersect {params.out}
		kmc_tools transform {params.out} dump {output.text}
		"""
# rule db_union:
# 	input:
# 		brain_kmc = base_path + "{study_id}/count_kmers/{genomes}_all_kmers.kmc_suf",
# 		genome_kmc = base_path + "{study_id}/count_kmers/31mers/{genomes}_all_kmers.kmc_suf"
# 	output:
# 		base_path + "{study_id}/count_kmers/{genomes}_union.kmc_suf"
# 	params:
# 		brain =base_path + "{study_id}/count_kmers/{genomes}_all_kmers",
# 		genome = base_path + "{study_id}/count_kmers/31mers/{genomes}_all_kmers",
# 		union = base_path + "{study_id}/count_kmers/{genomes}_union"
# 	log:
# 		base_path + "{study_id}/count_kmers/{genomes}_union_kmers.log"
# 	threads: 8
# 	resources:
# 		runtime = 240,
# 		mem_mb = 30000
# 	shell:
# 		"""
# 			kmc_tools simple {params.brain} {params.genome} union {params.union}
# 		"""

# rule db_subtract_mouse:
# 	input:
# 		brain_kmc = base_path + "{study_id}/count_kmers/mouse_all_kmers.kmc_suf",
# 		union = base_path + "{study_id}/count_kmers/human_union.kmc_suf"
# 	output:
# 		base_path + "{study_id}/databases/mouse_reset.kmc_suf"
# 	params:
# 		brain =base_path + "{study_id}/count_kmers/mouse_all_kmers",
# 		union = base_path + "{study_id}/count_kmers/human_union",
# 		subtract = base_path + "{study_id}/count_kmers/mouse_subtract",
# 		out =base_path + "{study_id}/databases/mouse_reset"
# 	log:
# 		base_path + "{study_id}/count_kmers/mouse_union_kmers.log"
# 	threads: 8
# 	resources:
# 		runtime = 240,
# 		mem_mb = 30000
# 	shell:
# 		"""
# 			kmc_tools simple {params.brain} {params.union} kmers_subtract {params.subtract}
# 			kmc_tools transform {params.subtract} set_counts 1 {params.out}
# 		"""

# rule db_subtract_human:
# 	input:
# 		brain_kmc = base_path + "{study_id}/count_kmers/human_all_kmers.kmc_suf",
# 		union = base_path + "{study_id}/count_kmers/mouse_union.kmc_suf"
# 	output:
# 		base_path + "{study_id}/databases/human_reset.kmc_suf"
# 	params:
# 		brain =base_path + "{study_id}/count_kmers/human_all_kmers",
# 		union = base_path + "{study_id}/count_kmers/mouse_union",
# 		subtract = base_path + "{study_id}/count_kmers/human_subtract",
# 		out =base_path + "{study_id}/databases/human_reset"
# 	log:
# 		base_path + "{study_id}/count_kmers/human_union_kmers.log"
# 	threads: 8
# 	resources:
# 		runtime = 240,
# 		mem_mb = 30000
# 	shell:
# 		"""
# 			kmc_tools simple {params.brain} {params.union} kmers_subtract {params.subtract}
# 			kmc_tools transform {params.subtract} set_counts 1 {params.out}
# 		"""

rule intersect_GBM43_kmers:
	input:
		GBM43gDNA = "/scratch/smankame/Murine/count_kmers/GBM43_gDNA_all_kmers.kmc_suf",
		GBMTMZgDNA= base_path + "{study_id}/count_kmers/human_unique_kmers.kmc_suf"
	output:
		base_path + "{study_id}/databases/GBM43_intersect_reset.kmc_suf"
	params:
		GBM43gDNA = "/scratch/smankame/Murine/count_kmers/GBM43_gDNA_all_kmers",
		GBMTMZgDNA = base_path + "{study_id}/count_kmers/human_unique_kmers",
		output_db = base_path + "{study_id}/databases/GBM43_intersect",
		reset = base_path + "{study_id}/databases/GBM43_intersect_reset"
	threads:
		16
	resources:
		runtime = 480,
		mem_mb = 50000
	shell:
		"""
			kmc_tools simple {params.GBM43gDNA} {params.GBMTMZgDNA} intersect {params.output_db}
			kmc_tools transform {params.output_db} set_counts 1 {params.reset}
		"""
# rule reset_dbs:
# 	input:
# 		kmc = base_path + "{study_id}/count_kmers/{genomes}_unique_kmers.kmc_suf",
# 	output:
# 		kmc = base_path + "{study_id}/count_kmers/{genomes}_unique_kmers_reset.kmc_suf",
# 	threads:
# 		16
# 	resources:
# 		runtime = 480,
# 		mem_mb = 50000
# 	params:
# 		inputs=base_path + "{study_id}/count_kmers/{genomes}_unique_kmers",
# 		outputs=base_path + "{study_id}/count_kmers/{genomes}_unique_kmers_reset"
# 	shell:
# 		"""
# 		kmc_tools transform {params.inputs} set_counts 1 {params.outputs}
# 		"""  

rule fastq_to_fasta:
	input:
		r1 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p1.fastq.gz",
		r2 = base_path + "{study_id}/bam_processing/ExtractUmis/{patient_id}/{sampleid}_processed.p2.fastq.gz"
	output:
		base_path + "{study_id}/fasta/{patient_id}/{sampleid}.fasta"
	resources:
		runtime = 480,
		mem_mb = 50000
	conda:
		"/home/smankame/miniforge3/envs/py2env"
	shell:
		"seqtk seq -A {input.r1} {input.r2} > {output}"

rule get_counts:
	input:
		fasta = base_path + "{study_id}/fasta/{patient_id}/{sampleid}.fasta",
		dbs = expand(base_path + "{study_id}/count_kmers/{genomes}_unique_kmers_reset.kmc_suf", study_id = study_filter, genomes = genomes)
	output:
		txt = (base_path + "{study_id}/get_counts/{patient_id}/{sampleid}.txt")
	log:
		base_path + "{study_id}/get_counts/{patient_id}/{sampleid}.log"
	params:
		db_mouse = base_path + "{study_id}/count_kmers/mouse_unique_kmers_reset",
		db_human = base_path + "{study_id}/databases/GBM43_intersect_reset"
	resources:
		runtime = 18500,
		mem_mb = 150000
	shell:
		"/tgen_labs/barthel/software/github/barthel/cfDNA/sWGS/workflow/scripts/get_counts_matrix {input.fasta} {params.db_human} {params.db_mouse} {output.txt} &>> {log}"


rule classify_reads_matrix:
	input:
		(base_path + "{study_id}/get_counts/{patient_id}/{sampleid}.txt")
	params:
		sample_id = "{sampleid}",
	output:
		base_path + "{study_id}/classified_reads/{patient_id}/{sampleid}.txt",
		base_path + "{study_id}/classified_reads/{patient_id}/{sampleid}_human.txt",
		base_path + "{study_id}/classified_reads/{patient_id}/{sampleid}_mouse.txt",
		base_path + "{study_id}/classified_reads/{patient_id}/{sampleid}_ambiguous.txt",
	threads:
		16
	resources:
		runtime = 480,
		mem_mb = 150000
	script:
		"/tgen_labs/barthel/software/github/barthel/cfDNA_sWGS_public/workflow/scripts/count_murine.py"
