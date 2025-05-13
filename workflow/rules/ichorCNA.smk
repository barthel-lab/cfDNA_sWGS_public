rule wig_files:
	input:
		bam = base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{sampleid}_BQSR_hg38.bam"
	output:
		wig =  base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.wig"
	params:
		index = base_path + "{study_id}/bam_processing/BQSR/{patient_id}/{sampleid}_BQSR_hg38.bam.bai",
		window_size = 1000000,
		quality = 20,
		chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
	shell:
		"""
		if [ ! -f {params.index} ]; then
			samtools index -b {input.bam}
		fi
		/home/smankame/miniforge3/envs/ichorcna/bin/readCounter --window {params.window_size} --quality {params.quality} \
			--chromosome {params.chromosomes} {input.bam} > {output.wig}
		"""

rule create_PON:
	input:
		pbmc = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC.wig"
	output:
		base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC_PON_median.rds" 
	params:
		file_list = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PON_filelist.txt",
		pon = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC_PON"
	shell:
		"""
		echo {input.pbmc} > {params.file_list}
		Rscript /home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/scripts/createPanelOfNormals.R \
			--filelist {params.file_list} \
			--gcWig /home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
			--centromere /home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
			--outfile {params.pon} 
		"""
rule run_ichorCNA:
	input:
		wig = base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.wig",
		PON = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC_PON_median.rds" 
		#PON = "/home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
	output:
		base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.cna.seg"
	params:
		outDir = str(base_path + "{study_id}/ichorcna/{patient_id}/"),
		gcWig = "/home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/inst/extdata/gc_hg38_1000kb.wig",
		mapWig = "/home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/inst/extdata/map_hg38_1000kb.wig",
		centro = "/home/smankame/miniforge3/envs/ichorcna/bin/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
		sample = "{sampleid}"
	script:
		"/tgen_labs/barthel/Brown_COH/scripts/ichorCNA/IchorCNA.sh"
