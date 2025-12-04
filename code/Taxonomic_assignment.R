#### SETTING THE ENVIRONMENT (and inputs)

suppressPackageStartupMessages({
  library(dada2)
  library(DECIPHER)
})

main <- function(input_paths = list(), output_paths = list(), params = list()) {
  script_path <- getwd()
  db <- params[["database"]]
  dada <- params[["dada_files_dir"]]
  out_dir <- params[["results_dir"]]

  db <- file.path(script_path, db)
  dada_path <- file.path(script_path, dada)
  output_path <- file.path(script_path, out_dir)

  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  setwd(dada_path)
  load(file = "seqtab.RData")

  ## TAXONOMY
  setwd(script_path)
  load(db)

  dna <- DNAStringSet(getSequences(seqtab.nochim))
  tax_info <- IdTaxa(test = dna, trainingSet = trainingSet, strand = "both", processors = NULL, threshold = 50)

  # RENAMING ASVs
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
  for (i in seq_len(dim(seqtab.nochim)[2])) {
    asv_headers[i] <- paste("ASV", i, sep = "_")
  }

  lengths <- vapply(seq_along(asv_seqs), function(i) nchar(asv_seqs[i]), integer(1))
  lengths_df <- data.frame(names = as.vector(asv_headers), lengths = as.vector(lengths))

  setwd(output_path)

  # WRITING FASTA
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  write(asv_fasta, "ASVs.fa")

  # WRITING ABUNDANCE TABLE
  asv_tab <- t(seqtab.nochim)
  row.names(asv_tab) <- sub(">", "", asv_headers)
  write.table(asv_tab, "ASVs_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)

  # WRITING TAX TABLE
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  asv_tax <- t(sapply(tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa
  }))
  colnames(asv_tax) <- ranks
  rownames(asv_tax) <- gsub(pattern = ">", replacement = "", x = asv_headers)

  write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)
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
  main(
    input_paths = list(),
    output_paths = list(),
    params = list(
      database = "db/SILVA_SSU_r138_2019.RData",
      dada_files_dir = ".",
      results_dir = "results_debug"
    )
  )
}