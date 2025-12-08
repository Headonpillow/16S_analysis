##### SETTING THE ENVIRONMENT (and inputs)

suppressPackageStartupMessages({
  library(dada2)
  library(gridExtra)
  library(dplyr)
})

main <- function(input_paths = list(), output_paths = list(), params = list()) {
  # Map Snakemake inputs
  samples_file <- input_paths[["samples"]]                    # intermediate/trimmed/samples.txt

  # Map Snakemake outputs
  track_out <- output_paths[["track"]]                        # results/denoising/read_count_tracking.tsv
  qc_out <- output_paths[["qc"]]                              # results/denoising/qc.pdf
  seqtab_out <- output_paths[["seqtab"]]                      # results/denoising/seqtab.RData

  # Map Snakemake params
  output_path <- params[["results_dir"]]                      # results/denoising 
  intermediate_path <- params[["intermediate_filtered_dir"]]  # intermediate/dada_filtered

  # Basic checks
  if (is.null(samples_file)) {
    stop("Missing required input: 'samples'")
  }
  if (any(vapply(list(track_out, qc_out, seqtab_out), is.null, logical(1)))) {
    stop("Missing required outputs: expect names 'track','qc','seqtab'.")
  }

  # Resolve paths
  # Use the samples directory as general input path (it should have been created from previous scripts)
  input_path <- dirname(samples_file)
  # Ensure the samples directory and file exist before proceeding
  if (!dir.exists(input_path)) {
    stop(sprintf("Samples directory not found: %s", input_path))
  }
  if (!file.exists(samples_file)) {
    stop(sprintf("Samples file not found: %s", samples_file))
  }

  # Read samples file
  samples <- scan(file = samples_file, what = "character")

  # Deduplicate sample names if needed (filterAndTrim requires distinct output filenames)
  if (any(duplicated(samples))) {
    warning(sprintf("Found %d duplicate sample name(s) in %s; keeping first occurrence of each sample.", sum(duplicated(samples)), samples_file))
    samples <- unique(samples)
  }

  # FILTERING
  # Build full paths for input fastq files (they live in the samples directory)
  forward_reads <- file.path(input_path, paste0(samples, ".R1.fastq.gz"))
  reverse_reads <- file.path(input_path, paste0(samples, ".R2.fastq.gz"))

  # Prepare filtered output file paths in the results directory (do not change cwd)
  # Ensure output directory exists
  out_dir <- if (!is.null(output_path) && nzchar(output_path)) output_path else dirname(seqtab_out)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  filtered_dir <- if (!is.null(intermediate_path) && nzchar(intermediate_path)) intermediate_path else out_dir
  dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
  filtered_forward_reads <- file.path(filtered_dir, paste0(samples, ".filtered.R1.fastq.gz"))
  filtered_reverse_reads <- file.path(filtered_dir, paste0(samples, ".filtered.R2.fastq.gz"))
  
  # QUALITY PLOTS
  quality_plot <- function(data) {
    plot <- plotQualityProfile(data)
    return(plot)
  }

  fplots <- lapply(forward_reads, quality_plot)
  rplots <- lapply(reverse_reads, quality_plot)

  # FILTERING
  filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, multithread = TRUE)

  filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]
  filtered_reverse_reads <- filtered_reverse_reads[file.exists(filtered_reverse_reads)]

  # DENOISE
  err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE)
  err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE)

  # TODO: pooling?
  dada_forward <- dada(filtered_forward_reads, err = err_forward_reads, pool = TRUE, multithread = TRUE)
  dada_reverse <- dada(filtered_reverse_reads, err = err_reverse_reads, pool = TRUE, multithread = TRUE)

  # MERGE
  merged_amplicons <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse, filtered_reverse_reads, verbose = TRUE)

  seqtab <- makeSequenceTable(merged_amplicons)

  # REMOVE SHORT SEQS
  seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 350:470]

  # CHIMERA REMOVAL
  seqtab.nochim <- removeBimeraDenovo(seqtab2, verbose = TRUE, method = "consensus", multithread = TRUE)

  ########### Outputs directory prepared above; do not change working directory
  # Ensure declared Snakemake output directories exist
  if (!is.null(track_out) && nzchar(track_out)) dir.create(dirname(track_out), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(qc_out) && nzchar(qc_out)) dir.create(dirname(qc_out), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(seqtab_out) && nzchar(seqtab_out)) dir.create(dirname(seqtab_out), recursive = TRUE, showWarnings = FALSE)

  # PERFORMANCE ASSESSMENT
  getN <- function(x) sum(getUniques(x))

  if (length(samples) != length(sapply(dada_reverse, getN))) {
    original <- file.path(out_dir, paste0(samples, ".filtered.R1.fastq.gz"))
    i <- match(original[!file.exists(original)], original)
    cat("sample: ", samples[i], " has been excluded from the analysis\n")
    samples <- samples[-i]
    filtered_out <- filtered_out[-i, ]
    rm(original, i)
  }

  summary_tab <- data.frame(row.names = samples,
                            dada2_input = filtered_out[, 1],
                            filtered = filtered_out[, 2],
                            dada_f = sapply(dada_forward, getN),
                            dada_r = sapply(dada_reverse, getN),
                            merged = sapply(merged_amplicons, getN),
                            nonchim = rowSums(seqtab.nochim),
                            final_perc_reads_retained = round(rowSums(seqtab.nochim)/filtered_out[, 1] * 100, 1))

  # Write the summary table to the declared Snakemake output path
  if (!is.null(track_out) && nzchar(track_out)) {
    write.table(summary_tab, track_out, quote = FALSE, sep = "\t", col.names = NA)
  } else {
    write.table(summary_tab, file.path(out_dir, "read_count_tracking.tsv"), quote = FALSE, sep = "\t", col.names = NA)
  }

  ################ OUTPUTS
  # SEQTAB
  if (!is.null(seqtab_out) && nzchar(seqtab_out)) {
    save(seqtab.nochim, file = seqtab_out)
  } else {
    save(seqtab.nochim, file = file.path(out_dir, "seqtab.RData"))
  }

  # GENERATE PDF OF QUALITY PLOTS
  len <- length(fplots)
  max_pages <- len/3
  if (len %% 3 != 0) {
    max_pages <- floor(max_pages)
  }

  i <- 1
  pdf(if (!is.null(qc_out) && nzchar(qc_out)) qc_out else file.path(out_dir, "qc.pdf"), onefile = TRUE, paper = "a4", height = 9, width = 9)
  while (i <= max_pages * 3) {
    z <- i + 2
    grid.arrange(grobs = c(fplots[i:z], rplots[i:z]), ncol = 2, as.table = FALSE)
    i <- i + 3
  }

  if (len %% 3 == 1) {
    grid.arrange(grobs = c(fplots[i], rplots[i]), ncol = 2, as.table = FALSE)
  } else if (len %% 3 == 2) {
    z <- i + 1
    grid.arrange(grobs = c(fplots[i:z], rplots[i:z]), ncol = 2, as.table = FALSE)
  }

  dev.off()
}

# --- Wrapper for Snakemake vs interactive execution ---
if (exists("snakemake")) {
  main(
    input_paths = as.list(snakemake@input),
    output_paths = as.list(snakemake@output),
    params = as.list(snakemake@params)
  )
} else {
    # Fallback for interactive debugging
    stop("This script is meant to be run by Snakemake. For interactive use, source() it and call main(...) manually.")
}
