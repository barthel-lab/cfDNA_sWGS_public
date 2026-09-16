rule wig_files:
	input:
		bam = bam_file_path + ".bam"
	output:
		wig =  base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.wig"
	params:
		index = bam_file_path + ".bam.bai",
		window_size = 1000000,
		quality = 20,
		chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY",
		path = readCounter_path
	resources:
		runtime = 200,
		mem_mb = 8000
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
	resources:
		runtime = 200,
		mem_mb = 8000
	shell:
		"""
		Rscript {params.ichor_path}/scripts/createPanelOfNormals.R \
			--filelist {params.file_list} \
			--gcWig {params.ichor_path}/inst/extdata/gc_hg38_1000kb.wig \
			--centromere {params.ichor_path}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
			--outfile {params.pon} 
		"""

##echo {input.pbmc} > {params.file_list}
rule run_ichorCNA:
	input:
		wig = base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.wig",
		PON = base_path + "{study_id}/ichorcna/Murine_PON_median.rds" 
		#PON = ichorcna_path+ "/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds" ###use if patient normal is not available
	output:
		base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.cna.seg"
	params:
		outDir = str(base_path + "{study_id}/ichorcna/{patient_id}/"),
		gcWig = ichorcna_path+"/inst/extdata/gc_hg38_1000kb.wig",
		mapWig = ichorcna_path+"/inst/extdata/map_hg38_1000kb.wig",
		centro = ichorcna_path+"/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
		sample = "{sampleid}",
		ichor_path = ichorcna_path
	log:
		base_path + "{study_id}/ichorcna/{patient_id}/{sampleid}.log"
	resources:
		runtime = 200,
		mem_mb = 10000
	shell:
		"""
			Rscript {params.ichor_path}/scripts/runIchorCNA.R \
				--id {params.sample} \
				--WIG {input.wig} \
				--ploidy "c(2,3)" \
				--normal "c(0.95, 0.99, 0.995, 0.999)" \
				--maxCN 3 \
				--gcWig {params.gcWig} \
				--mapWig {params.mapWig} \
				--centromere {params.centro} \
				--normalPanel {input.PON} \
				--includeHOMD False \
				--chrs "c(1:22)" \
				--chrTrain "c(1:22)" \
				--estimateNormal True \
				--estimatePloidy False \
				--estimateScPrevalence False \
				--scStates "c()" \
				--txnE 0.9999999 \
				--txnStrength 100000000 \
				--outDir {params.outDir} &> {log}
		"""

def get_patient_seg_files(wildcards):
    """Get all seg files for a given patient within a study."""
    patient_df = sWGS_table[
        (sWGS_table['Study'] == wildcards.study_id) & 
        (sWGS_table['Patient'] == wildcards.patient_id)
    ]
    samples = patient_df['Sample'].tolist()
    return expand(
        file_path + "/{study_id}/ichorcna/{patient_id}/{sampleid}.cna.seg",
        study_id=wildcards.study_id,
        patient_id=wildcards.patient_id,
        sampleid=samples
    )

rule arm_level_aneuploidy_per_patient:
    input:
        seg = get_patient_seg_files
    output:
        summary = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/taylor_aneuploidy_summary.txt",
        taylor  = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/taylor_aneuploidy.txt",
        prop    = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/prop_aneuploidy.txt",
        cnv_stats  = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/CNV_stats.txt",
        seg_stats  = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/seg_stats.txt",
        heatmap    = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/heatmap.pdf"
    params:
        output_dir = file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}"
    log:
        file_path + "/{study_id}/arm_level_aneuploidy/{patient_id}/aneuploidy.log"
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "/tgen_labs/barthel/software/github/barthel/cfDNA_sWGS_public/workflow/scripts/aneuploidy_ichorcna.R"

#### GATK pipeline (I've found this to overcall for cfDNA, it's more optimized for tumor samples)
# rule addReadGroups:
# 	input:
# 		bam = bam_file_path
# 	output:
# 		bam = base_path + "{study_id}/gatk_cnv/readGroups/{patient_id}/{sampleid}_RG_hg38.bam"
# 	log:
# 		base_path + "{study_id}/gatk_cnv/readGroups/{patient_id}/{sampleid}_RG.log"
# 	message:
# 		"Adding RG information to sample {wildcards.sampleid}."
# 	params:
# 		paramID = "{sampleid}"
# 		picard_fp = picard_path
# 	threads: 8
# 	resources:
# 		runtime = 200,
# 		mem_mb = 15000
# 	shell:
# 		"""
# 		java -Xmx45g -jar {params.picard_fp} AddOrReplaceReadGroups \
# 			--INPUT {input.bam} \
# 			--OUTPUT {output.bam} \
# 			--RGLB {params.paramID} \
# 			--RGPL ILLUMINA \
# 			--RGPU {params.paramID} \
# 			--RGSM {params.paramID}
# 		"""

# rule sort_index_bam:
#     input:
#         base_path + "{study_id}/gatk_cnv/readGroups/{patient_id}/{sampleid}_RG_hg38.bam"
#     output:
#         bam=base_path + "{study_id}/gatk_cnv/readGroups/{patient_id}/{sampleid}_RG_sorted.bam",
#         bai=base_path + "{study_id}/gatk_cnv/readGroups/{patient_id}/{sampleid}_RG_sorted.bam.bai"
#     shell:
#         """
#         samtools sort -o {output.bam} {input}
#         samtools index {output.bam}
#         """
# rule preprocessIntervals:
# 	input:
# 		reference = human
# 	output:
# 		base_path +"references/GRCh38/Homo_sapiens_assembly38.100kbp.interval_list"
# 	shell:
# 		"""
# 			gatk PreprocessIntervals \
# 				-R {input.reference} \
# 				--bin-length 100000 \
# 				--padding 0 \
# 				-O {output}
# 		"""

# rule collectreadcounts:
# 	input:
# 		bam = bam_fp,
# 		interv = gatk_intervals_100kbp ##change based on which interval size user wants
# 	output:
# 		base_path + "{study_id}/gatk_cnv/CollectReadCounts/{patient_id}/{sampleid}_counts.hdf5"
# 	threads: 8
# 	resources:
# 		mem_mb = 50000
# 	shell:
# 		"""
# 		gatk --java-options -Xmx8g CollectReadCounts \
# 			-I {input.bam} \
# 			-L {input.interv} \
# 			--interval-merging-rule OVERLAPPING_ONLY \
# 			-O {output}
# 		"""

# # Create a pon
# # Comment this if you want to use reference genome as background of copy number
# # rule CreateReadCountPanelOfNormals:
# # 	input:
# # 		hdf5 = base_path + "{study_id}/gatk_cnv/CollectReadCounts/{patient_id}/{patient_id}_PBMC_counts.hdf5"
# # 	output:
# # 		base_path + "{study_id}/gatk_cnv/CreateReadCountPanelOfNormals/{patient_id}.pon.hdf5"
# # 	threads: 8
# # 	resources:
# # 		mem_mb = 50000
# # 	shell:
# # 		"""
# # 		gatk --java-options -Xmx12g CreateReadCountPanelOfNormals \
# # 			-I {input.hdf5} \
# # 			-O {output}
# # 		"""

# # denoise reads counts
# # Comment input.pon if you want to use reference genome as background of copy number
# # Disable --count-panel-of-normals option if you want to use reference genome as background of copy number
# rule DenoiseReadCounts:
# 	input:
# 		hdf5 = base_path + "{study_id}/gatk_cnv/CollectReadCounts/{patient_id}/{sampleid}_counts.hdf5",
# 		#pon = base_path + "{study_id}/gatk_cnv/CreateReadCountPanelOfNormals/{patient_id}.pon.hdf5 --count-panel-of-normals {input.pon} \
# 	output:
# 		standardized = base_path + "{study_id}/gatk_cnv/DenoiseReadCounts/{patient_id}/{sampleid}.standardizedCR.tsv",
# 		denoised = base_path + "{study_id}/gatk_cnv/DenoiseReadCounts/{patient_id}/{sampleid}.denoisedCR.tsv"
# 	shell:
# 		"""
# 		gatk --java-options -Xmx12g DenoiseReadCounts \
# 			-I {input.hdf5} \
# 			--standardized-copy-ratios {output.standardized} \
# 			--denoised-copy-ratios {output.denoised}
# 		"""

# # plot denoised and standardized CNVsq
# rule PlotDenoisedCopyRatios:
# 	input:
# 		standardized =base_path + "{study_id}/gatk_cnv/DenoiseReadCounts/{patient_id}/{sampleid}.standardizedCR.tsv",
# 		denoised = base_path + "{study_id}/gatk_cnv/DenoiseReadCounts/{patient_id}/{sampleid}.denoisedCR.tsv",
# 		gatk_dict = hg38_dict
# 	output:
# 		base_path + "{study_id}/gatk_cnv/plotcr/{patient_id}/{sampleid}.denoised.png",
# 		base_path + "{study_id}/gatk_cnv/plotcr/{patient_id}/{sampleid}.standardizedMAD.txt",
# 		base_path + "{study_id}/gatk_cnv/plotcr/{patient_id}/{sampleid}.denoisedMAD.txt",
# 		base_path + "{study_id}/gatk_cnv/plotcr/{patient_id}/{sampleid}.deltaMAD.txt",
# 		base_path + "{study_id}/gatk_cnv/plotcr/{patient_id}/{sampleid}.scaledDeltaMAD.txt"
# 	params:
# 		outputdir = base_path + "{study_id}/gatk_cnv/plotcr/{patient_id}/{sampleid}",
# 		outputprefix = "{sampleid}",
# 		minChrSize = "46709982"
# 	conda:
# 		"/home/ychen/miniconda3/envs/R_v3.6"
# 	shell:
# 		"""
# 		gatk --java-options -Xmx12g PlotDenoisedCopyRatios \
# 			--standardized-copy-ratios {input.standardized} \
# 			--denoised-copy-ratios {input.denoised} \
# 			--sequence-dictionary {input.gatk_dict} \
# 			--minimum-contig-length {params.minChrSize} \
# 			--output {params.outputdir} \
# 			--output-prefix {params.outputprefix}
# 		"""

# rule CollectAllelicCounts:
# 	input:
# 		bam =  bam_fp,
# 		fasta = human,
# 		vcf = vcf_intervals
# 	output:
# 	   base_path + "{study_id}/gatk_cnv/CollectAllelicCounts/{patient_id}/{sampleid}.allelicCounts.tsv",
# 	resources:
# 		mem_gb=80
# 	shell:
# 		"""
# 		gatk --java-options "-Xmx60g" CollectAllelicCounts \
# 			-I {input.bam} \
# 			-R {input.fasta} \
# 			-L {input.vcf} \
# 			-O {output}
# 		"""

# rule ModelSegments:
# 	input:
# 		denoised = base_path + "{study_id}/gatk_cnv/DenoiseReadCounts/{patient_id}/{sampleid}.denoisedCR.tsv",
# 		allelicCounts = base_path + "{study_id}/gatk_cnv/CollectAllelicCounts/{patient_id}/{sampleid}.allelicCounts.tsv"
# 	output:
# 		base_path + "{study_id}/gatk_cnv/ModelSegments/{patient_id}/{sampleid}.cr.seg"
# 	params:
# 		outdir = base_path + "{study_id}/gatk_cnv/ModelSegments/{patient_id}/",
# 		outprefix = "{sampleid}"
# 	threads: 20
# 	resources:
# 		mem_mb=300000
# 	shell:
# 		"""
# 		gatk --java-options "-Xmx300g" ModelSegments \
# 			--denoised-copy-ratios {input.denoised} \
# 			--allelic-counts {input.allelicCounts} \
# 			--output-prefix {params.outprefix} \
# 			--output {params.outdir}
# 		"""

# rule CallCopyRatioSegments:
# 	input:
# 		base_path + "{study_id}/gatk_cnv/ModelSegments/{patient_id}/{sampleid}.cr.seg",
# 	output:
# 		base_path + "{study_id}/gatk_cnv/ModelSegments/{patient_id}/{sampleid}.called.seg"
# 	resources:
# 		mem_mb=80000
# 	shell:
# 		"""
# 		gatk --java-options "-Xmx80g"  CallCopyRatioSegments \
# 			--input {input} \
# 			--output {output}
# 		"""

# # BUG from the gatk, it is not working with CallCopyRatioSegments
# # rule PlotModeledSegments:
# #     input:
# #         denoised = "results/DenoiseReadCounts/{sample}/{sample}.denoisedCR.tsv",
# #         allelicCounts = "results/CollectAllelicCounts/{sample}/hg38_{sample}.allelicCounts.tsv",
# #         segments = "results/ModelSegments/{sample}/{sample}.cr.seg",
# #         gatk_dict = dictionary
# #     output:
# #         "results/PlotModeledSegments/{sample}/{sample}_PlotModeledSegments.png"
# #     params:
# #         outprefix = "{sample}",
# #         outdir = "results/PlotModeledSegments/{sample}",
# #         minChromSize = "46709982"
# #     conda:
# #         "/home/ychen/miniconda3/envs/R_v3.6_dup" # Strange bug, R_v3.6_dup and R_v3.6 are the same env, but R_v3.6_dup works.
# #     shell:"""
# #         gatk --java-options "-Xmx150g" PlotModeledSegments \
# #             --denoised-copy-ratios {input.denoised} \
# #             --allelic-counts {input.allelicCounts} \
# #             --segments {input.segments} \
# #             --sequence-dictionary {input.gatk_dict} \
# #             --minimum-contig-length {params.minChromSize} \
# #             --output {params.outdir} \
# #             --output-prefix {params.outprefix}
# #         """
# # need to figure out why some column is missing in the input