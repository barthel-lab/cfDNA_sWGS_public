library(XenofilteR)

# --- Snakemake inputs/params ---
human_bam <- snakemake@input[["human_bam"]]
mouse_bam  <- snakemake@input[["mouse_bam"]]
sampleid   <- snakemake@wildcards[["sampleid"]]
outdir     <- snakemake@params[["outdir"]]
n_workers  <- snakemake@threads

# Create directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

sample.list <- data.frame(
  graft = human_bam,
  host  = mouse_bam,
  stringsAsFactors = FALSE
)

bp.param <- SnowParam(
  workers = n_workers,
  type = "SOCK"
)

XenofilteR(
  sample.list        = sample.list,
  destination.folder = outdir,
  bp.param           = bp.param,
  output.names       = sampleid
)