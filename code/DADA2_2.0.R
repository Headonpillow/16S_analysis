#####SETTING THE ENVIRONMENT (and inputs)

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))

# Setting up the different input and output paths, retrieving them from Snakemake
# once taken the directories from snakemake params, the full path is created by pasting them
# on the path of the current directory (from where the Snakefile is executed)
script_path <- getwd()
in_dir <- snakemake@params[["sample_file_loc"]]
out_dir <- snakemake@params[["results_dir"]]
input_path <- paste0(script_path, "/", in_dir)
output_path <- paste0(script_path, "/", out_dir)

# Setting the input directory to the correct path
setwd(input_path)

samples <- scan(file="samples.txt", what="character")

#FILTERING

forward_reads <- paste0(samples, ".R1.fastq.gz")
reverse_reads <- paste0(samples, ".R2.fastq.gz")

filtered_forward_reads <- paste0(samples, ".filtered.R1.fastq.gz")
filtered_reverse_reads <- paste0(samples, ".filtered.R2.fastq.gz")

#QUALITY PLOTS

quality_plot <- function(data) {
  plot <- plotQualityProfile(data)
  return(plot)
}

fplots <- lapply(forward_reads, quality_plot)
rplots <- lapply(reverse_reads, quality_plot)

#FILTERING

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, multithread=TRUE)

filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]
filtered_reverse_reads <- filtered_reverse_reads[file.exists(filtered_reverse_reads)]

#DENOISE

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

#plotErrors(err_forward_reads, nominalQ=TRUE)
#plotErrors(err_reverse_reads, nominalQ=TRUE)

# TODO: pooling?
dada_forward <- dada(filtered_forward_reads, err=err_forward_reads, pool=TRUE, multithread=TRUE)
dada_reverse <- dada(filtered_reverse_reads, err=err_reverse_reads, pool=TRUE, multithread=TRUE)

#MERGE

merged_amplicons <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse, filtered_reverse_reads, verbose=TRUE)

seqtab <- makeSequenceTable(merged_amplicons)

#REMOVE SHORT SEQS
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 410:470]

#CHIMERA REMOVAL

seqtab.nochim <- removeBimeraDenovo(seqtab2, verbose=T, method = "consensus", multithread = TRUE)

########### Set working directory to output path
setwd(output_path)

#PERFORMANCE ASSESSMENT

getN <- function(x) sum(getUniques(x))

if(length(samples) != length(sapply(dada_reverse, getN))) {
  original <- paste0(samples, ".filtered.R1.fastq.gz")
  i <- match(original[!file.exists(original)], original)
  cat("sample: ", samples[i], " has been excluded from the analysis")
  samples <- samples[-i]
  filtered_out <- filtered_out[-i,]
  rm (original, i)
}

summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1], filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN), dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN), nonchim=rowSums(seqtab.nochim), final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100,1))

write.table(summary_tab, "read_count_tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

################OUTPUTS

#SEQTAB

save(seqtab.nochim, file="./seqtab.RData")

#GENERATE PDF OF QUALITY PLOTS

len <- length(fplots)
max_pages <- len/3

if(len%%3 != 0) {
  max_pages <- floor(max_pages)
}

i <- 1

pdf("qc.pdf", onefile = TRUE, paper = "a4", height = 9, width = 9)

while (i <= max_pages*3) {
  z <- i+2
  grid.arrange( grobs = c(fplots[i: z], rplots[i: z]), ncol = 2, as.table = FALSE)
  i <- i+3
}

if(len%%3 == 1) {
  grid.arrange( grobs = c(fplots[i], rplots[i]), ncol = 2, as.table = FALSE)
} else if (len%%3 == 2) {
  z <- i+1
  grid.arrange( grobs = c(fplots[i: z], rplots[i: z]), ncol = 2, as.table = FALSE)
}

dev.off()
