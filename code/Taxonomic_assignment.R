#### SETTING THE ENVIRONMENT (and inputs)

suppressPackageStartupMessages({
  library(dada2)
  library(DECIPHER)
})

main <- function(input_paths = list(), output_paths = list(), params = list()) {
  # Map Snakemake inputs
  seqtab_file <- input_paths[["seqtab"]]    # results/denoising/seqtab.RData
  samples_file <- input_paths[["samples"]]  # intermediate/trimmed/samples.txt

  # Map Snakemake outputs
  fa_out <- output_paths[["fa"]]            # results/asv/ASVs.fa
  counts_out <- output_paths[["counts"]]    # results/asv/ASVs_counts.tsv
  tax_out <- output_paths[["tax"]]          # results/asv/ASVs_taxonomy.tsv

  # Map Snakemake params
  db_param <- params[["database"]]

  # Use the directory of the fa_out as the general output dir
  output_dir <- dirname(fa_out)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Basic checks
  if (is.null(seqtab_file)) {
    stop("Missing required input: 'seqtab'")
  }
  if (any(vapply(list(is.null(fa_out), is.null(counts_out), is.null(tax_out)), isTRUE, logical(1)))) {
    stop("Missing required outputs: expect names 'fa','counts','tax'.")
  }

  # Resolve path for database
  db_path <- db_param
  # Ensure paths are valid
  if (!is.null(db_path) && !file.exists(db_path)) {
    db_path <- file.path(getwd(), db_param)
  }

  # Load seqtab and database (do not change working directory)
  load(seqtab_file)
  if (!is.null(db_path) && file.exists(db_path)) load(db_path)

  # TAXONOMY
  dna <- DNAStringSet(getSequences(seqtab.nochim))
  tax_info <- IdTaxa(test = dna, trainingSet = trainingSet, strand = "both", processors = NULL, threshold = 50)

  # RENAMING ASVs
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
  for (i in seq_len(dim(seqtab.nochim)[2])) {
    asv_headers[i] <- paste(">ASV", i, sep = "_")
  }

  lengths <- vapply(seq_along(asv_seqs), function(i) nchar(asv_seqs[i]), integer(1))
  lengths_df <- data.frame(names = as.vector(asv_headers), lengths = as.vector(lengths))

  # WRITING FASTA
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  write(asv_fasta, fa_out)

  # WRITING ABUNDANCE TABLE
  asv_tab <- t(seqtab.nochim)
  row.names(asv_tab) <- sub(">", "", asv_headers)
  write.table(asv_tab, counts_out, sep = "\t", quote = FALSE, col.names = NA)

  # WRITING TAX TABLE
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  asv_tax <- t(sapply(tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa
  }))
  colnames(asv_tax) <- ranks
  rownames(asv_tax) <- gsub(pattern = ">", replacement = "", x = asv_headers)

  write.table(asv_tax, tax_out, sep = "\t", quote = FALSE, col.names = NA)
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
