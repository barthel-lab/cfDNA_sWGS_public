rule multisample_mutect2:
	input:
		normal = base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{patient_id}_PBMC_BQSR_hg38.bam",
		germ = germline,
		fasta = hg38,
		pon = nPON,
		vcf = hotspot_file
	params:
		patient = "{patient_id}",
		sample = lambda wildcards: " ".join(
			f"-I {base_path}{s_id}/bam_processing/BQSR/{p_id}/{sample_id}_BQSR_hg38.bam"
			for s_id, p_id, sample_id in zip(study_list, patient_list, sample_list)
			if p_id == wildcards.patient_id and s_id == wildcards.study_id
		)
	output:
		vcf = base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample.vcf.gz",
		f1r2=base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample_f1r2.tar.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	shell:
		"""
			gatk Mutect2 \
				-R {input.fasta} \
				{params.sample} \
				-alleles {input.vcf} \
				-normal {params.patient}_PBMC \
				--germline-resource {input.germ} \
				--genotype-germline-sites true --genotype-pon-sites true \
				--f1r2-tar-gz {output.f1r2} \
				--panel-of-normals {input.pon} \
				-O {output.vcf}
		"""

###-alleles {input.vcf} \

rule LearnReadOrientationModel:
	input:
		base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample_f1r2.tar.gz"
	output:
		base_path + "{study_id}/LearnReadOrientationModel/{patient_id}/{patient_id}_multisample_read-orientation-model.tar.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/LearnReadOrientationModel/{patient_id}.LearnReadOrientationModel.log"
	shell:
		"""
			gatk LearnReadOrientationModel \
				-I {input} -O {output} \
				> {log} 2>&1
	"""
#run on every sample
rule GetPileupSummaries:
	input:
		bam = base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{sampleid}_BQSR_hg38.bam",
		var = variant,
		inter = intervals,
		fasta = hg38
	output:
		base_path + "{study_id}/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries_table.txt"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/GetPileupSummaries/{patient_id}/{sampleid}.GetPileupSummaries.log"
	shell:
		"""
			gatk GetPileupSummaries -I {input.bam} -V {input.var} \
				-L {input.inter}  -O {output} \
				-R {input.fasta} > {log} 2>&1
		"""
##run only for tumor sample with matched normal mode
rule CalculateContamination:
	input:
		tumor_pile =base_path + "{study_id}/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries_table.txt",
		normal_pile = base_path + "{study_id}/GetPileupSummaries/{patient_id}/{patient_id}_PBMC_GetPileupSummaries_table.txt"
	output:
		contam = base_path + "{study_id}/CalculateContamination/{patient_id}/{sampleid}_contamination.txt",
		seg = base_path + "{study_id}/CalculateContamination/{patient_id}/{sampleid}_segmentation.tsv"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/GetPileupSummaries/{patient_id}/{sampleid}_CalculateContamination.log"
	shell:
		"""
			gatk CalculateContamination -I {input.tumor_pile} \
				-matched {input.normal_pile} \
				--tumor-segmentation {output.seg} \
				-O {output.contam} > {log} 2>&1
		"""
				
rule FilterMutectCalls:
	input:
		vcf=base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample.vcf.gz",
		fasta = hg38
	params:
		contam = lambda wildcards: " ".join(f"--contamination-table {base_path}/{s_id}/CalculateContamination/{p_id}/{sample_id}_contamination.txt" for s_id, p_id, sample_id in zip(study_list, patient_list, sample_list)
			if p_id == wildcards.patient_id and s_id == wildcards.study_id),
		seg = lambda wildcards: " ".join(f"--tumor-segmentation {base_path}/{s_id}/CalculateContamination/{p_id}/{sample_id}_segmentation.tsv" for s_id, p_id, sample_id in zip(study_list, patient_list, sample_list)
			if p_id == wildcards.patient_id and s_id == wildcards.study_id)
	output:
		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_filtered.vcf.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/mutect2/filtered_{patient_id}.var.mutect2.log"
	shell:
		"""
			gatk FilterMutectCalls -R {input.fasta} -V {input.vcf} \
				{params.contam} \
				{params.seg} \
				--max-events-in-region 10 \
				--f-score-beta 1.5  \
				-O {output} > {log} 2>&1
		"""

rule removeFailedCalls:
	input: 
		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_filtered.vcf.gz"
	output:
		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_pass.vcf.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/mutect2/pass_{patient_id}.var.mutect2.log"
	shell:
		"""
			(bcftools view -f "PASS" {input} -Oz -o {output}) > {log} 2>&1
			tabix -f -p vcf {output}
		"""
rule funcotator:
	input:
		vcf = base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_pass.vcf.gz",
		fasta = hg38,
		FUN = fun_lib
	output:
		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_anno.vcf.gz"
	params:
		genome = "hg38",
		FORMAT = "VCF"
	shell:
		"""
			gatk Funcotator \
				--variant {input.vcf} \
				--reference {input.fasta} \
				--ref-version {params.genome} \
				--data-sources-path {input.FUN} \
				--output {output} \
				--output-file-format {params.FORMAT} 
		"""

