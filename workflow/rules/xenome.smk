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


