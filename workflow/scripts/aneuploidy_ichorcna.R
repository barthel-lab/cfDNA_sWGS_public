#!/usr/bin/env Rscript
# ============================================================================
# Arm-Level Aneuploidy Analysis (Taylor Aneuploidy Score) - ichorCNA Version
# ============================================================================
# v2 (corrected) — adapted from arm_level_aneuploidy_ichorcna.R
#
# CHANGE vs original (the "all sub-threshold / score 0" fix):
#   The original called an arm gain/loss only if the SINGLE LONGEST CONTIGUOUS
#   run of same-direction segments covered >80% of the arm (prop_arm built from
#   `partion_id == 1`). A whole-arm event broken into pieces by one noisy bin was
#   split into short runs, none reaching 80%, and fell through to NA. Highly
#   aneuploid-but-fragmented genomes (e.g. UPN210, prop_aneuploidy 0.68) scored 0.
#
#   FIX: measure prop_arm as the TOTAL fraction of the arm covered by each call
#   direction (summed across all runs), then call the arm by its dominant
#   direction. The >0.80 arm-level threshold is UNCHANGED — only the quantity it
#   is applied to is corrected. Set ARM_CALL_THRESHOLD below to tune.
#
#   Also self-contained: arm assignment uses a base-R centromere split instead of
#   GenomicRanges, and patchwork (dendrogram) is optional.
#
# Usage:
#   Rscript arm_level_aneuploidy_ichorcna_v2.R <input(s) or dir> [output_dir]
# ============================================================================

ARM_CALL_THRESHOLD <- 0.80   # fraction of the arm that must be one direction

# ============================================================================
# HARDCODED ARM LENGTHS (hg38)
# ============================================================================
get_arm_lengths <- function() {
  arm_data <- data.frame(
    chrom = c(
      "chr1", "chr1", "chr2", "chr2", "chr3", "chr3", "chr4", "chr4",
      "chr5", "chr5", "chr6", "chr6", "chr7", "chr7", "chr8", "chr8",
      "chr9", "chr9", "chr10", "chr10", "chr11", "chr11", "chr12", "chr12",
      "chr13", "chr13", "chr14", "chr14", "chr15", "chr15", "chr16", "chr16",
      "chr17", "chr17", "chr18", "chr18", "chr19", "chr19", "chr20", "chr20",
      "chr21", "chr21", "chr22", "chr22", "chrX", "chrX", "chrY", "chrY"
    ),
    arm = c(
      "p", "q", "p", "q", "p", "q", "p", "q",
      "p", "q", "p", "q", "p", "q", "p", "q",
      "p", "q", "p", "q", "p", "q", "p", "q",
      "p", "q", "p", "q", "p", "q", "p", "q",
      "p", "q", "p", "q", "p", "q", "p", "q",
      "p", "q", "p", "q", "p", "q", "p", "q"
    ),
    start = c(
      0, 123400000, 0, 93900000, 0, 90900000, 0, 50000000,
      0, 48800000, 0, 59800000, 0, 60100000, 0, 45200000,
      0, 43000000, 0, 39800000, 0, 53400000, 0, 35500000,
      0, 17700000, 0, 17200000, 0, 19000000, 0, 36800000,
      0, 25100000, 0, 18500000, 0, 26200000, 0, 28100000,
      0, 12000000, 0, 15000000, 0, 61000000, 0, 10400000
    ),
    end = c(
      123400000, 248956422, 93900000, 242193529, 90900000, 198295559, 50000000, 190214555,
      48800000, 181538259, 59800000, 170805979, 60100000, 159345973, 45200000, 145138636,
      43000000, 138394717, 39800000, 133797422, 53400000, 135086622, 35500000, 133275309,
      17700000, 114364328, 17200000, 107043718, 19000000, 101991189, 36800000, 90338345,
      25100000, 83257441, 18500000, 80373285, 26200000, 58617616, 28100000, 64444167,
      12000000, 46709983, 15000000, 50818468, 61000000, 156040895, 10400000, 57227415
    ),
    stringsAsFactors = FALSE
  )
  arm_data$length <- arm_data$end - arm_data$start
  return(arm_data)
}

# ============================================================================
# PACKAGE CHECK  (GenomicRanges no longer required; patchwork optional)
# ============================================================================
check_packages <- function() {
  required_packages <- c("dplyr", "purrr", "tidyr", "ggplot2", "stringr")
  missing <- c()
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) missing <- c(missing, pkg)
  }
  if (length(missing) > 0) {
    message("ERROR: Missing required R packages:")
    for (pkg in missing) message(sprintf("  - %s (install.packages('%s'))", pkg, pkg))
    stop("Please install missing packages and re-run.", call. = FALSE)
  }
  message("All required packages available.")
}
check_packages()

suppressMessages({
  library(dplyr); library(purrr); library(tidyr)
  library(ggplot2); library(stringr)
})
HAVE_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

# ============================================================================
# INPUT PARSING
# ============================================================================
extract_sample_id <- function(filepath) {
  sample_id <- basename(filepath)
  sample_id <- sub("\\.cna\\.seg$", "", sample_id)
  sample_id <- sub("\\.seg$", "", sample_id)
  return(sample_id)
}
find_seg_files <- function(dir_path) {
  list.files(dir_path, pattern = "\\.(cna\\.)?seg$", full.names = TRUE)
}

if (exists("snakemake")) {
  message("Running in Snakemake mode")
  seg_files <- snakemake@input[["seg"]]
  if (is.null(seg_files)) seg_files <- unlist(snakemake@input)
  seg_files <- as.character(seg_files)
  output_dir <- snakemake@params[["output_dir"]]
  sample_ids <- sapply(seg_files, extract_sample_id)
  multi_sample <- length(seg_files) > 1
  output_prefix <- if (multi_sample) "combined" else sample_ids[1]
} else {
  message("Running in CLI mode")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript ...R <input> [output_dir]", call. = FALSE)
  seg_files <- c(); output_dir <- getwd()
  for (i in seq_along(args)) {
    arg <- args[i]
    is_last <- (i == length(args))
    is_seg_file <- grepl("\\.(cna\\.)?seg$", arg, ignore.case = TRUE)
    is_existing_dir <- dir.exists(arg)
    if (is_existing_dir && !is_seg_file) {
      if (i == 1) {
        found_files <- find_seg_files(arg)
        if (length(found_files) > 0) seg_files <- c(seg_files, found_files)
        else stop(sprintf("No .seg files found in directory: %s", arg), call. = FALSE)
      } else output_dir <- arg
    } else if (is_seg_file && file.exists(arg)) {
      seg_files <- c(seg_files, arg)
    } else if (is_last && !is_seg_file && length(seg_files) > 0) {
      output_dir <- arg
    } else if (!file.exists(arg) && is_seg_file) {
      stop(sprintf("File not found: %s", arg), call. = FALSE)
    } else if (!is_seg_file && !is_existing_dir) {
      output_dir <- arg
    }
  }
  if (length(seg_files) == 0) stop("No valid seg files found.", call. = FALSE)
  sample_ids <- sapply(seg_files, extract_sample_id)
  multi_sample <- length(seg_files) > 1
  output_prefix <- if (multi_sample) "combined" else sample_ids[1]
}

message(sprintf("Number of samples: %d", length(seg_files)))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# READ ICHORCNA SEGMENTATION FILE
# ============================================================================
read_ichorcna_seg <- function(seg_file, sample_id) {
  df <- read.delim(seg_file, header = TRUE, stringsAsFactors = FALSE)
  logR_col <- grep("\\.logR$", colnames(df), value = TRUE)
  if (length(logR_col) == 0) {
    if ("median" %in% colnames(df)) logR_col <- "median"
    else stop(sprintf("No logR/median column in %s", seg_file), call. = FALSE)
  }
  colnames(df)[colnames(df) == "chr"]   <- "CONTIG"
  colnames(df)[colnames(df) == "start"] <- "START"
  colnames(df)[colnames(df) == "end"]   <- "END"
  colnames(df)[colnames(df) == logR_col] <- "MEAN_LOG2_COPY_RATIO"
  df$SampleID <- sample_id
  if (!any(grepl("^chr", df$CONTIG))) df$CONTIG <- paste0("chr", df$CONTIG)
  df
}

seg_list <- list()
for (i in seq_along(seg_files)) {
  seg_list[[i]] <- read_ichorcna_seg(seg_files[i], sample_ids[i])
  message(sprintf("Loaded %d segments for %s", nrow(seg_list[[i]]), sample_ids[i]))
}
gatk_seg <- bind_rows(seg_list)
arm_lengths <- get_arm_lengths()

# ============================================================================
# NEUTRAL-SEGMENT WEIGHTED STATS  (unchanged from original)
# ============================================================================
unfiltered_seg_wmean_wsd <- gatk_seg %>%
  mutate(COPY_RATIO = 2^MEAN_LOG2_COPY_RATIO, SEG_LENGTH = END - START + 1) %>%
  filter(COPY_RATIO >= 0.9 & COPY_RATIO <= 1.1) %>%
  group_by(SampleID) %>%
  summarise(num_seg = n(),
            wmean = sum(SEG_LENGTH * COPY_RATIO) / sum(SEG_LENGTH),
            wsd = ifelse(num_seg > 1,
                         sqrt(sum(SEG_LENGTH * (COPY_RATIO - wmean)^2) / sum(SEG_LENGTH)),
                         NA_real_),
            .groups = "drop")

gatk_seg <- gatk_seg %>%
  mutate(COPY_RATIO = 2^MEAN_LOG2_COPY_RATIO, SEG_LENGTH = END - START + 1)
gatk_seg_joined <- gatk_seg %>% left_join(unfiltered_seg_wmean_wsd, by = "SampleID")

filtered_seg_wmean_wsd <- gatk_seg_joined %>%
  filter(COPY_RATIO >= 0.9 & COPY_RATIO <= 1.1,
         abs(COPY_RATIO - wmean) <= 2 * wsd | is.na(wsd)) %>%
  group_by(SampleID) %>%
  summarise(fnum_seg = n(),
            fwmean = sum(SEG_LENGTH * COPY_RATIO) / sum(SEG_LENGTH),
            fwsd = ifelse(
              fnum_seg > 1 & sum(SEG_LENGTH * (COPY_RATIO^2)) > (sum(SEG_LENGTH * COPY_RATIO)^2 / sum(SEG_LENGTH)),
              sqrt((sum(SEG_LENGTH * COPY_RATIO^2) - (sum(SEG_LENGTH * COPY_RATIO)^2) / sum(SEG_LENGTH)) /
                     (sum(SEG_LENGTH) - 1)),
              0),
            .groups = "drop")

# ============================================================================
# SEGMENT GAIN/LOSS/NEUTRAL CALL  (unchanged from original)
# ============================================================================
gatk_seg_with_stats <- gatk_seg %>% left_join(filtered_seg_wmean_wsd, by = "SampleID")
gatk_seg_with_cnv_call <- gatk_seg_with_stats %>%
  mutate(CNV_CALL = case_when(
    COPY_RATIO >= 0.9 & COPY_RATIO <= 1.1 ~ 0L,
    (COPY_RATIO - fwmean) < -2.0 * fwsd ~ -1L,
    (COPY_RATIO - fwmean) >  2.0 * fwsd ~  1L,
    TRUE ~ 0L))

# ---- loss segment stats (unchanged) ----
unfiltered_loss_wmean_wsd <- gatk_seg_with_cnv_call %>% filter(CNV_CALL == -1) %>%
  group_by(SampleID) %>%
  summarise(num_seg = n(),
            wmean = sum((END - START + 1) * COPY_RATIO) / sum(END - START + 1),
            wsd = sqrt((sum((END - START + 1) * COPY_RATIO^2) -
                          (sum((END - START + 1) * COPY_RATIO)^2 / sum(END - START + 1))) /
                         (sum(END - START + 1) - 1)),
            .groups = 'drop') %>%
  mutate(wsd = ifelse(is.na(wsd) | wsd < 0, NA_real_, wsd))
if (nrow(unfiltered_loss_wmean_wsd) > 0) {
  filtered_loss_wmean_wsd <- gatk_seg_with_cnv_call %>%
    inner_join(unfiltered_loss_wmean_wsd, by = "SampleID") %>%
    filter(CNV_CALL == -1, abs(COPY_RATIO - wmean) <= 2 * wsd | is.na(wsd)) %>%
    group_by(SampleID) %>%
    summarise(fnum_seg = n(),
              fwmean = sum((END - START + 1) * COPY_RATIO) / sum(END - START + 1),
              fwsd = ifelse(fnum_seg > 1,
                            sqrt(pmax(sum((END - START + 1) * (COPY_RATIO^2)) -
                                        (sum((END - START + 1) * COPY_RATIO)^2) / sum(END - START + 1)) /
                                   (sum(END - START + 1) - 1)), NA_real_),
              .groups = "drop") %>%
    mutate(fwsd = ifelse(is.na(fwsd) | is.infinite(fwsd), NA_real_, fwsd))
} else {
  filtered_loss_wmean_wsd <- data.frame(SampleID=character(0), fnum_seg=integer(0),
                                        fwmean=numeric(0), fwsd=numeric(0))
}
# ---- gain segment stats (unchanged) ----
unfiltered_gain_wmean_wsd <- gatk_seg_with_cnv_call %>% filter(CNV_CALL == 1) %>%
  group_by(SampleID) %>%
  summarise(num_seg = n(),
            wmean = sum((END - START - 1) * COPY_RATIO) / sum(END - START - 1),
            wsd = ifelse(num_seg > 1,
                         sqrt(pmax((sum((END - START - 1) * (COPY_RATIO^2)) -
                                      (sum((END - START - 1) * COPY_RATIO)^2) / sum(END - START - 1)) /
                                     (sum(END - START - 1) - 1), 0)), NA_real_),
            .groups = 'drop') %>%
  mutate(wsd = ifelse(is.na(wsd) | is.infinite(wsd), NA_real_, wsd))
if (nrow(unfiltered_gain_wmean_wsd) > 0) {
  filtered_gain_wmean_wsd <- gatk_seg_with_cnv_call %>%
    inner_join(unfiltered_gain_wmean_wsd, by = "SampleID") %>%
    filter(CNV_CALL == 1, abs(COPY_RATIO - wmean) <= 2 * wsd | is.na(wsd)) %>%
    group_by(SampleID) %>%
    summarise(fnum_seg = n(),
              fwmean = sum((END - START - 1) * COPY_RATIO) / sum(END - START - 1),
              fwsd = ifelse(fnum_seg > 1,
                            sqrt(pmax(sum((END - START - 1) * (COPY_RATIO^2)) -
                                        (sum((END - START - 1) * COPY_RATIO)^2) / sum(END - START - 1)) /
                                   (sum(END - START - 1) - 1)), NA_real_),
              .groups = "drop") %>%
    mutate(fwsd = ifelse(is.na(fwsd) | is.infinite(fwsd), NA_real_, fwsd))
} else {
  filtered_gain_wmean_wsd <- data.frame(SampleID=character(0), fnum_seg=integer(0),
                                        fwmean=numeric(0), fwsd=numeric(0))
}

seg_stats <- gatk_seg_with_cnv_call %>% select(SampleID) %>% distinct() %>%
  left_join(filtered_seg_wmean_wsd, by = "SampleID") %>%
  left_join(filtered_loss_wmean_wsd, by = "SampleID", suffix = c("", "_loss")) %>%
  left_join(filtered_gain_wmean_wsd, by = "SampleID", suffix = c("", "_gain"))

# ============================================================================
# ARM ASSIGNMENT  (base-R centromere split; replaces GenomicRanges)
# ============================================================================
# Each chromosome has one centromere = end of its p arm. A segment is on the p
# arm if it starts before the centromere, else the q arm. This matches the
# original findOverlaps + !duplicated (first overlap = p) behaviour.
centromere <- arm_lengths %>% filter(arm == "p") %>% transmute(seqnames = chrom, cen = end)

seg_df <- gatk_seg_with_cnv_call %>%
  transmute(seqnames = CONTIG, start = START, end = END,
            CNV_CALL, COPY_RATIO = 2^MEAN_LOG2_COPY_RATIO, SampleID) %>%
  left_join(centromere, by = "seqnames") %>%
  mutate(arm = ifelse(start < cen, "p", "q")) %>%
  select(-cen)
seg_df$seg_size <- seg_df$end - seg_df$start + 1
seg_df <- merge(seg_df, arm_lengths[, c("chrom", "arm", "length")],
                by.x = c("seqnames", "arm"), by.y = c("chrom", "arm"), all.x = TRUE)
colnames(seg_df)[colnames(seg_df) == "length"] <- "arm_size"

# ============================================================================
# ARM-LEVEL CALL  (FIX: total per-direction coverage, not longest contiguous run)
# ============================================================================
# prop_arm denominator is the CALLABLE arm length (sum of binned/segmented sizes
# in the arm), NOT the full hg38 cytogenetic arm length. ichorCNA only bins the
# mappable genome, so 15-25% of each cytogenetic arm (centromeres, telomeres,
# low-mappability) is never measured; dividing by the full arm length caps a
# fully-altered arm below 100% and can make a near-complete event unreachable
# (e.g. chr21q is only ~75% mappable -> can never hit 80% of the full arm).
# Dividing by callable length asks "of what was measured, how much is altered".
message(sprintf("Calling arms by per-direction coverage of CALLABLE arm (threshold %.2f)...",
                ARM_CALL_THRESHOLD))

# Callable arm length = total size of segments with a usable copy ratio in the arm.
callable_arm <- seg_df %>%
  filter(COPY_RATIO > 0, seg_size > 0) %>%
  group_by(SampleID, seqnames, arm) %>%
  summarise(callable_size = sum(seg_size),
            arm_size = min(arm_size, na.rm = TRUE),
            .groups = "drop")

arm_direction_coverage <- seg_df %>%
  filter(COPY_RATIO > 0, seg_size > 0) %>%
  group_by(SampleID, seqnames, arm, CNV_CALL) %>%
  summarise(dir_size = sum(seg_size),
            num_seg  = n(),
            wcr      = sum(seg_size * COPY_RATIO) / sum(seg_size),
            msq      = sum(seg_size * COPY_RATIO^2) / sum(seg_size),
            .groups  = "drop") %>%
  left_join(callable_arm, by = c("SampleID", "seqnames", "arm")) %>%
  mutate(wsd = sqrt(pmax(msq - wcr^2, 0)),
         prop_arm = dir_size / callable_size,            # <- callable denominator
         callable_frac = callable_size / arm_size)

cnv_by_arm <- arm_direction_coverage %>%
  group_by(SampleID, seqnames, arm) %>%
  arrange(desc(prop_arm), .by_group = TRUE) %>%
  slice(1) %>%                       # dominant direction by total arm coverage
  ungroup() %>%
  mutate(arm_call = case_when(
           CNV_CALL ==  1 & prop_arm > ARM_CALL_THRESHOLD ~  1,
           CNV_CALL == -1 & prop_arm > ARM_CALL_THRESHOLD ~ -1,
           CNV_CALL ==  0                                 ~  0,
           TRUE ~ NA_real_),
         arm_num_seg  = num_seg,
         arm_cr_wmean = wcr,
         arm_cr_wsd   = wsd) %>%
  select(SampleID, seqnames, arm, arm_call, arm_num_seg, arm_cr_wmean, arm_cr_wsd)

# ============================================================================
# TAYLOR ANEUPLOIDY SCORES + PROPORTION  (unchanged from original)
# ============================================================================
taylor_aneuploidy <- cnv_by_arm %>%
  group_by(SampleID) %>%
  summarise(
    aneuploidy_score = sum(!is.na(arm_call) & arm_call != 0),
    aneuploidy_amp_score = sum(!is.na(arm_call) & arm_call == 1),
    aneuploidy_del_score = sum(!is.na(arm_call) & arm_call == -1),
    max_loss_arm_n = ifelse(all(is.na(arm_cr_wmean)), NA, arm_num_seg[which.min(arm_cr_wmean)]),
    max_loss_arm_wmean = min(arm_cr_wmean, na.rm = TRUE),
    max_loss_arm_wsd = ifelse(all(is.na(arm_cr_wmean)), NA, arm_cr_wsd[which.min(arm_cr_wmean)]),
    max_gain_arm_n = ifelse(all(is.na(arm_cr_wmean)), NA, arm_num_seg[which.max(arm_cr_wmean)]),
    max_gain_arm_wmean = max(arm_cr_wmean, na.rm = TRUE),
    max_gain_arm_wsd = ifelse(all(is.na(arm_cr_wmean)), NA, arm_cr_wsd[which.max(arm_cr_wmean)]),
    .groups = "drop")

prop_aneuploidy <- gatk_seg_with_cnv_call %>%
  filter(!CONTIG %in% c("chrX", "chrY")) %>%
  group_by(SampleID) %>%
  summarise(seg_size = sum(END - START - 1),
            het_size = sum(ifelse(CNV_CALL == 0, END - START - 1, 0)),
            .groups = "drop")

taylor_aneuploidy_summary <- taylor_aneuploidy %>%
  left_join(prop_aneuploidy, by = "SampleID") %>%
  mutate(prop_aneuploidy = round(1.0 - het_size / seg_size, 4)) %>%
  select(SampleID, prop_aneuploidy, aneuploidy_score, aneuploidy_amp_score,
         aneuploidy_del_score, max_loss_arm_n, max_loss_arm_wmean, max_loss_arm_wsd,
         max_gain_arm_n, max_gain_arm_wmean, max_gain_arm_wsd)

# ============================================================================
# SAVE OUTPUT FILES
# ============================================================================
message("Saving output files...")
write.table(prop_aneuploidy, file.path(output_dir, paste0(output_prefix, "_prop_aneuploidy.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(taylor_aneuploidy_summary, file.path(output_dir, paste0(output_prefix, "_taylor_aneuploidy_summary.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(taylor_aneuploidy, file.path(output_dir, paste0(output_prefix, "_taylor_aneuploidy.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(cnv_by_arm, file.path(output_dir, paste0(output_prefix, "_CNV_stats.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(seg_stats, file.path(output_dir, paste0(output_prefix, "_seg_stats.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# ============================================================================
# HEATMAP  (simple tile map; dendrogram only if patchwork is available)
# ============================================================================
standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
chr_arm_seg_data <- cnv_by_arm %>%
  filter(seqnames %in% standard_chroms, !is.na(arm)) %>%
  mutate(chr_arm = paste0(seqnames, arm)) %>%
  select(SampleID, chr_arm, arm_call)

if (nrow(chr_arm_seg_data) > 0) {
  heatmap_long <- chr_arm_seg_data %>%
    mutate(arm_call_label = case_when(
      arm_call == -1 ~ "Loss", arm_call == 0 ~ "Neutral", arm_call == 1 ~ "Gain",
      TRUE ~ "Sub-threshold (<80%)")) %>%
    mutate(arm_call_label = factor(arm_call_label,
                                   levels = c("Loss", "Neutral", "Gain", "Sub-threshold (<80%)")))
  arm_order <- paste0("chr", rep(1:22, each = 2), c("p", "q"))
  arm_order <- arm_order[arm_order %in% unique(heatmap_long$chr_arm)]
  heatmap_long$chr_arm <- factor(heatmap_long$chr_arm, levels = arm_order)

  heatmap_plot <- ggplot(heatmap_long, aes(x = SampleID, y = chr_arm, fill = arm_call_label)) +
    geom_tile(color = "grey70", width = 1, height = 1) +
    scale_fill_manual(values = c("Loss"="blue","Neutral"="white","Gain"="red",
                                 "Sub-threshold (<80%)"="grey80"),
                      name = "Arm Call", drop = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
          axis.text.y = element_text(size = 4), axis.title = element_blank(),
          panel.grid = element_blank()) +
    labs(title = sprintf("Arm-Level Aneuploidy - ichorCNA v2 (%d samples)",
                         length(unique(heatmap_long$SampleID))))
  n_samples <- length(unique(heatmap_long$SampleID)); n_arms <- length(arm_order)
  pdf(file.path(output_dir, paste0(output_prefix, "_heatmap.pdf")),
      width = max(1, n_samples * 0.25 + 2), height = max(1, n_arms * 0.25 + 2))
  suppressWarnings(print(heatmap_plot)); dev.off()
}

message(sprintf("Done. Output prefix '%s' in %s", output_prefix, output_dir))
