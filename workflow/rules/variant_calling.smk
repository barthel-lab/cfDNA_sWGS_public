# rule multisample_mutect2:
# 	input:
# 		normal = base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{patient_id}_PBMC_BQSR_hg38.bam",
# 		germ = germline,
# 		fasta = hg38,
# 		pon = nPON,
# 		vcf = hotspot_file
# 	params:
# 		patient = "{patient_id}",
# 		sample = lambda wildcards: " ".join(
# 			f"-I {base_path}{s_id}/bam_processing/BQSR/{p_id}/{sample_id}_BQSR_hg38.bam"
# 			for s_id, p_id, sample_id in zip(study_list, patient_list, sample_list)
# 			if p_id == wildcards.patient_id and s_id == wildcards.study_id
# 		)
# 	output:
# 		vcf = base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample.vcf.gz",
# 		f1r2=base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample_f1r2.tar.gz"
# 	threads:
# 		16
# 	resources:
# 		mem_mb = 250000
# 	shell:
# 		"""
# 			gatk Mutect2 \
# 				-R {input.fasta} \
# 				{params.sample} \
# 				-alleles {input.vcf} \
# 				-normal {params.patient}_PBMC \
# 				--germline-resource {input.germ} \
# 				--genotype-germline-sites true --genotype-pon-sites true \
# 				--f1r2-tar-gz {output.f1r2} \
# 				--panel-of-normals {input.pon} \
# 				-O {output.vcf}
# 		"""

# ###-alleles {input.vcf} \

# rule LearnReadOrientationModel:
# 	input:
# 		base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample_f1r2.tar.gz"
# 	output:
# 		base_path + "{study_id}/LearnReadOrientationModel/{patient_id}/{patient_id}_multisample_read-orientation-model.tar.gz"
# 	threads:
# 		16
# 	resources:
# 		mem_mb = 250000
# 	log:
# 		base_path + "{study_id}/logs/LearnReadOrientationModel/{patient_id}.LearnReadOrientationModel.log"
# 	shell:
# 		"""
# 			gatk LearnReadOrientationModel \
# 				-I {input} -O {output} \
# 				> {log} 2>&1
# 	"""
# #run on every sample
# rule GetPileupSummaries:
# 	input:
# 		bam = base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{sampleid}_BQSR_hg38.bam",
# 		var = variant,
# 		inter = intervals,
# 		fasta = hg38
# 	output:
# 		base_path + "{study_id}/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries_table.txt"
# 	threads:
# 		16
# 	resources:
# 		mem_mb = 250000
# 	log:
# 		base_path + "{study_id}/logs/GetPileupSummaries/{patient_id}/{sampleid}.GetPileupSummaries.log"
# 	shell:
# 		"""
# 			gatk GetPileupSummaries -I {input.bam} -V {input.var} \
# 				-L {input.inter}  -O {output} \
# 				-R {input.fasta} > {log} 2>&1
# 		"""
# ##run only for tumor sample with matched normal mode
# rule CalculateContamination:
# 	input:
# 		tumor_pile =base_path + "{study_id}/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries_table.txt",
# 		normal_pile = base_path + "{study_id}/GetPileupSummaries/{patient_id}/{patient_id}_PBMC_GetPileupSummaries_table.txt"
# 	output:
# 		contam = base_path + "{study_id}/CalculateContamination/{patient_id}/{sampleid}_contamination.txt",
# 		seg = base_path + "{study_id}/CalculateContamination/{patient_id}/{sampleid}_segmentation.tsv"
# 	threads:
# 		16
# 	resources:
# 		mem_mb = 250000
# 	log:
# 		base_path + "{study_id}/logs/GetPileupSummaries/{patient_id}/{sampleid}_CalculateContamination.log"
# 	shell:
# 		"""
# 			gatk CalculateContamination -I {input.tumor_pile} \
# 				-matched {input.normal_pile} \
# 				--tumor-segmentation {output.seg} \
# 				-O {output.contam} > {log} 2>&1
# 		"""
				
# rule FilterMutectCalls:
# 	input:
# 		vcf=base_path + "{study_id}/mutect2/{patient_id}/{patient_id}_multisample.vcf.gz",
# 		fasta = hg38
# 	params:
# 		contam = lambda wildcards: " ".join(f"--contamination-table {base_path}/{s_id}/CalculateContamination/{p_id}/{sample_id}_contamination.txt" for s_id, p_id, sample_id in zip(study_list, patient_list, sample_list)
# 			if p_id == wildcards.patient_id and s_id == wildcards.study_id),
# 		seg = lambda wildcards: " ".join(f"--tumor-segmentation {base_path}/{s_id}/CalculateContamination/{p_id}/{sample_id}_segmentation.tsv" for s_id, p_id, sample_id in zip(study_list, patient_list, sample_list)
# 			if p_id == wildcards.patient_id and s_id == wildcards.study_id)
# 	output:
# 		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_filtered.vcf.gz"
# 	threads:
# 		16
# 	resources:
# 		mem_mb = 250000
# 	log:
# 		base_path + "{study_id}/logs/mutect2/filtered_{patient_id}.var.mutect2.log"
# 	shell:
# 		"""
# 			gatk FilterMutectCalls -R {input.fasta} -V {input.vcf} \
# 				{params.contam} \
# 				{params.seg} \
# 				--max-events-in-region 10 \
# 				--f-score-beta 1.5  \
# 				-O {output} > {log} 2>&1
# 		"""

# rule removeFailedCalls:
# 	input: 
# 		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_filtered.vcf.gz"
# 	output:
# 		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_pass.vcf.gz"
# 	threads:
# 		16
# 	resources:
# 		mem_mb = 250000
# 	log:
# 		base_path + "{study_id}/mutect2/pass_{patient_id}.var.mutect2.log"
# 	shell:
# 		"""
# 			(bcftools view -f "PASS" {input} -Oz -o {output}) > {log} 2>&1
# 			tabix -f -p vcf {output}
# 		"""
# rule funcotator:
# 	input:
# 		vcf = base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_pass.vcf.gz",
# 		fasta = hg38,
# 		FUN = fun_lib
# 	output:
# 		base_path +  "{study_id}/mutect2/{patient_id}/{patient_id}_anno.vcf.gz"
# 	params:
# 		genome = "hg38",
# 		FORMAT = "VCF"
# 	shell:
# 		"""
# 			gatk Funcotator \
# 				--variant {input.vcf} \
# 				--reference {input.fasta} \
# 				--ref-version {params.genome} \
# 				--data-sources-path {input.FUN} \
# 				--output {output} \
# 				--output-file-format {params.FORMAT} 
# 		"""

################### SINGLE SAMPLE
rule wgs_mutect2:
   input:
	   bam = base_path + "{study_id}/xenofilteR/{patient_id}/{sampleid}/Filtered_bams/{sampleid}_Filtered.bam",
	   fasta= human,
	   germ = germline,
	   pon = nPON,
	   inter = intervals,
	   vcf = "/tgen_labs/barthel/references/GRCh38/lifted_hotspot_idh1_tert_grc37.vcf"
   output:
	   unfVar =base_path +  "{study_id}/mutect2/{patient_id}/{sampleid}.vcf.gz",
	   orient = base_path + "{study_id}/mutect2/{patient_id}/{sampleid}_f1r2.tar.gz"
   threads:
	   16
   resources:
	   mem_mb = 250000
   log:
	   base_path + "{study_id}/logs/mutect2/{patient_id}/{sampleid}_var.mutect2.log"
   shell:
	   """
		   gatk Mutect2  \
		   -R {input.fasta} \
		   -I {input.bam} \
		   --f1r2-tar-gz {output.orient} \
		   --germline-resource {input.germ} \
		   --alleles {input.vcf} \
		   --genotype-germline-sites true --genotype-pon-sites true \
		   --panel-of-normals {input.pon} \
		   -O {output.unfVar} > {log} 2>&1
	   """
#
rule wgs_LearnReadOrientationModel:
	input:
		base_path + "{study_id}/mutect2/{patient_id}/{sampleid}_f1r2.tar.gz",
	output:
		base_path + "{study_id}/LearnReadOrientationModel/{patient_id}/{sampleid}_read-orientation-model.tar.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/LearnReadOrientationModel/{patient_id}/{sampleid}_LearnReadOrientationModel.log"
	shell:
		"""
			gatk LearnReadOrientationModel \
				-I {input} -O {output} \
				> {log} 2>&1
		"""

rule wgs_GetPileupSummaries:
	input:
		bam = base_path + "{study_id}/xenofilteR/{patient_id}/{sampleid}/Filtered_bams/{sampleid}_Filtered.bam",
		var = variant,
		inter = intervals,
		fasta = human
	output:
		base_path + "{study_id}/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries_table.txt"
	threads:
			16
	resources:
			mem_mb = 250000
	log:
			base_path + "{study_id}/logs/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries.log"
	shell:
		"""
			gatk GetPileupSummaries -I {input.bam} -V {input.var} \
				-L {input.inter} \
				-O {output} \
				-R {input.fasta} > {log} 2>&1
		"""
rule wgs_CalculateContamination:
	input:
		base_path + "{study_id}/GetPileupSummaries/{patient_id}/{sampleid}_GetPileupSummaries_table.txt"
	output:
		base_path + "{study_id}/CalculateContamination/{patient_id}/{sampleid}_contamination.txt"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/GetPileupSummaries/{patient_id}/{sampleid}_CalculateContamination.log"
	shell:
		"""
			gatk CalculateContamination -I {input} \
				-O {output} > {log} 2>&1
		"""

rule wgs_FilterMutectCalls:
	input:
		vcf=base_path +  "{study_id}/mutect2/{patient_id}/{sampleid}.vcf.gz",
		fasta = human
	params:
		contam=base_path + "{study_id}/CalculateContamination/{patient_id}/{sampleid}_contamination.txt",
		orientation =base_path + "{study_id}/LearnReadOrientationModel/{patient_id}/{sampleid}_read-orientation-model.tar.gz"
	output:
		base_path +  "{study_id}/mutect2/{patient_id}/{sampleid}_filtered.vcf.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/mutect2/{patient_id}/{sampleid}_filtered.log"
	shell:
		"""
			gatk FilterMutectCalls -R {input.fasta} -V {input.vcf} \
				--contamination-table {params.contam} \
				--ob-priors {params.orientation} \
				--max-events-in-region 5 \
				--f-score-beta 1.5  \
				-O {output} > {log} 2>&1
		"""

rule wgs_removeFailedCalls:
	input: 
		base_path +  "{study_id}/mutect2/{patient_id}/{sampleid}_filtered.vcf.gz"
	output:
		base_path + "{study_id}/mutect2/{patient_id}/{sampleid}_pass.vcf.gz"
	threads:
		16
	resources:
		mem_mb = 250000
	log:
		base_path + "{study_id}/logs/mutect2/{patient_id}/{sampleid}_pass.log"
	shell:
		"""
			(bcftools view -f "PASS" {input} -Oz -o {output}) > {log} 2>&1
			tabix -f -p vcf {output}
		"""

rule wgs_funcotator:
	input:
		vcf = base_path + "{study_id}/mutect2/{patient_id}/{sampleid}_pass.vcf.gz",
		fasta = human,
		FUN = fun_lib
	output:
		base_path + "{study_id}/mutect2/{patient_id}/{sampleid}_anno.vcf.gz"
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

rule wgs_extract_VAF:
	input:
		base_path + "{study_id}/mutect2/{patient_id}/{sampleid}_anno.vcf.gz"
	output:
		txt = base_path +  "{study_id}/mutect2/{patient_id}/{sampleid}_anno.txt"
	params:
		sample_list = lambda wildcards: wildcards.sampleid	
	shell:
		"""
			echo -e "CHROM\tPOS\tSAMPLE\tREF\tALT\tAF\tDP\tGENE\tTYPE\tVARIANT\tCHANGE" > {output.txt}

			bcftools query -s {params.sample_list} -f '%CHROM\t%POS\t[%AD]\t[%AF]\t[%DP]\t%INFO/FUNCOTATION\n' {input} | \
			awk -F '\t' '{{ 
				split($3, ad, ","); 
				split($6, a, "|"); 
				gsub(/\[/, "", a[1]);
				print $1"\t"$2"\t"{params.sample_list}"\t"ad[1]"\t"ad[2]"\t"$4"\t"$5"\t"a[1]"\t"a[6]"\t"a[8]"\t"a[19] 
			}}' >> {output.txt}	
		"""