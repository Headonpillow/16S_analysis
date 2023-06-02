####SETTING THE ENVIRONMENT (and inputs)

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(DECIPHER))

script_path <- getwd()
db <- snakemake@params[["database"]]
dada <- snakemake@params[["dada_files_dir"]]
out_dir <- snakemake@params[["results_dir"]]

db <- paste0(script_path, "/", db)
dada_path <- paste0(script_path, "/", dada)
output_path <- paste0(script_path, "/", out_dir)

setwd(dada_path)
load(file="seqtab.RData")

##TAXONOMY

setwd(script_path)
load(db)

dna <- DNAStringSet(getSequences(seqtab.nochim))
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="top", processors=NULL, threshold=60)

#RENAMING ASVs

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

lengths <- vector()
for (i in 1:length(asv_seqs)) {
  asv_len <- length(strsplit(asv_seqs, "")[[i]])
  lengths <- append(lengths, asv_len)
}
lengths_df <- data.frame(names = as.vector(asv_headers))
lengths_df <- cbind(lengths_df, lengths = as.vector(lengths))

setwd(output_path)

#WRITING FASTA

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#WRITING ABUNDANCE TABLE

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#WRITING TAX TABLE

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)