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

#FILTERING

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE=c(2,4), trimLeft=c(20,20), truncLen=c(290,220), rm.phix = TRUE, multithread=TRUE)

filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]
filtered_reverse_reads <- filtered_reverse_reads[file.exists(filtered_reverse_reads)]

#LEARN ERRORS

err_forward_reads <- learnErrors(filtered_forward_reads, nbases= 1e8, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, nbases= 1e8, multithread=TRUE)

#DENOISE AND MERGE

mergers <- vector("list", length(samples))
names(mergers) <- samples
dadas_f <- vector("list", length(samples))
names(dadas_f) <- samples
dadas_r <- vector("list", length(samples))
names(dadas_r) <- samples
for(sam in samples) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(paste0("./", sam, ".filtered.R1.fastq.gz"))
    ddF <- dada(derepF, err=err_forward_reads, multithread=TRUE)
    dadas_f[[sam]] <- ddF
    derepR <- derepFastq(paste0("./", sam, ".filtered.R2.fastq.gz"))
    ddR <- dada(derepR, err=err_reverse_reads, multithread=TRUE)
    dadas_r[[sam]] <- ddR
    merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)
    mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
seqtab <- makeSequenceTable(mergers)

#REMOVE SHORT SEQS
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 350:460]

#CHIMERA REMOVAL

seqtab.nochim <- removeBimeraDenovo(seqtab2, verbose=T, method = "consensus", multithread = TRUE)

########### Set working directory to output path
setwd(output_path)

#PERFORMANCE ASSESSMENT

getN <- function(x) sum(getUniques(x))

if(length(samples) != length(sapply(dadas_r, getN))) {
  original <- paste0(samples, ".filtered.R1.fastq.gz")
  i <- match(original[!file.exists(original)], original)
  cat("sample: ", samples[i], " has been excluded from the analysis")
  samples <- samples[-i]
  filtered_out <- filtered_out[-i,]
  rm (original, i)
}

summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1], filtered=filtered_out[,2], dada_f=sapply(dadas_f, getN), dada_r=sapply(dadas_r, getN), merged=sapply(mergers, getN), nonchim=rowSums(seqtab.nochim), final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100,1))

write.table(summary_tab, "read_count_tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

################OUTPUTS

#SEQTAB

save(seqtab.nochim, file="./seqtab.RData")
