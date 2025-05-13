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
		path = readCounter_path
	shell:
		"""
		if [ ! -f {params.index} ]; then
			samtools index -b {input.bam}
		fi
		{params.path} --window {params.window_size} --quality {params.quality} \
			--chromosome {params.chromosomes} {input.bam} > {output.wig}
		"""

rule create_PON:
	input:
		pbmc = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC.wig"
	output:
		base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC_PON_median.rds" 
	params:
		file_list = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PON_filelist.txt",
		pon = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC_PON",
		ichor_path = ichorcna_path
	shell:
		"""
		echo {input.pbmc} > {params.file_list}
		Rscript {params.ichor_path}/scripts/createPanelOfNormals.R \
			--filelist {params.file_list} \
			--gcWig {params.ichor_path}/inst/extdata/gc_hg38_1000kb.wig \
			--centromere {params.ichor_path}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
			--outfile {params.pon} 
		"""
rule run_ichorCNA:
	input:
		wig = base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.wig",
		PON = base_path + "{study_id}/ichorcna/{patient_id}/{patient_id}_PBMC_PON_median.rds" 
		#PON = ichorcna_path+"/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
	output:
		base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.cna.seg"
	params:
		outDir = str(base_path + "{study_id}/ichorcna/{patient_id}/"),
		gcWig = ichorcna_path+"/inst/extdata/gc_hg38_1000kb.wig",
		mapWig = ichorcna_path+"/inst/extdata/map_hg38_1000kb.wig",
		centro = ichorcna_path+"/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
		sample = "{sampleid}",
		ichor_path = ichorcna_path
	script:
    shell:
        """
			Rscript {params.ichor_path}/scripts/runIchorCNA.R \
				--id {params.sample} \
				--WIG {input.wig} \
				--ploidy "c(2,3)" \
				--normal "c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)" \
				--maxCN 3 \
				--gcWig {params.gcWig} \
				--mapWig {params.mapWig} \
				--centromere {params.centro} \
				--normalPanel {input.PON} \
				--includeHOMD False \
				--chrs "c(1:22)" \
				--chrTrain "c(1:22)" \
				--estimateNormal True \
				--estimatePloidy True \
				--estimateScPrevalence False \
				--scStates "c(1)" \
				--txnE 0.9999999 \
				--txnStrength 100000000 \
				--outDir {params.outDir}
        """
