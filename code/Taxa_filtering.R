############### ENVIRONMENT SETUP (loading libraries and all the snakemake params)

# Libraries
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(vegan))

# Snakemake paths and params
script_path <- getwd()
samples <- snakemake@params[["sample_file_loc"]]
asv_files <- snakemake@params[["asv_dir"]]
metadata <- snakemake@params[["metadata_dir"]]
out_dir <- snakemake@params[["results_dir"]]

samples_path <- paste0(script_path, "/", samples)
asv_files_path <- paste0(script_path, "/", asv_files)
metadata_path <- paste0(script_path, "/", metadata)
output_path <- paste0(script_path, "/", out_dir)

############### PHYLOSEQ OBJECT CREATION

# Start with loading sample names
setwd(samples_path)
samples <- scan(file="samples.txt", what="character")

# Now loading the ASVs count_table and the ASVs tax_table
setwd(asv_files_path)
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t", fill=TRUE)
# Replace colnames of count_table with the actual sample names
colnames(count_tab) <- samples
tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T, row.names=1, check.names=F, sep="\t"), fill=TRUE)
# Also, read the ASVs sequences for further processing and for selecting the NAs to blast
ASVs <- read.fasta("ASVs.fa")

# Load the metadata and check which ones match with the sample names of the files (keep only those)
setwd(metadata_path)
sample_info_tab <- read.table("metadata.tsv", header=T, row.names=1, check.names=F, sep="\t", fill=TRUE)
sample_info_tab <- filter(sample_info_tab, rownames(sample_info_tab) %in% samples)
order <- match(rownames(sample_info_tab), colnames(count_tab))
sample_info_tab <- arrange(sample_info_tab, order)

# Finally put all the informations together in a Phyloseq object
Phyloseq_object <- phyloseq(otu_table(count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab))

# Filter all the samples with 0 counts
zero_counts <- colSums(otu_table(Phyloseq_object))==0
Phyloseq_object <- prune_samples(!zero_counts, Phyloseq_object)
condition <- all(zero_counts == FALSE)

if(condition == TRUE){
    cat("\n")
    cat("No samples were removed.")
    cat("\n")
} else{
    cat("\n")
    filtered <- which(zero_counts == TRUE)
    cat(paste("WARN: sample", names(zero_counts[filtered]), "have been filtered out because it had 0 counts."))
    cat("\n")
}

############### DESEQ2

# Retrieve data from the filtered Phyloseq in data.frame format:
sample_info_tab <- data.frame(sample_data(Phyloseq_object))
count_tab <- data.frame(otu_table(Phyloseq_object))

# The formula has been set to "formula ~1" because this equals to no formula, since we are not using this data for differential expression.
metadata_fields <- colnames(data.frame(sample_data(Phyloseq_object)))
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~1)
# Since getting an error with the dataset mosquitoes+water:
# "Error in estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means"
# It's because the table is quite sparse with many zeroes, so I needed to add the
# Next command to deal with it.
# If there is not such problem, next line can be commented.
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)

# Now create a phyloseq object with normalized counts
Phyloseq_object_vst <- phyloseq(otu_table(vst_trans_count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab))

############### STARTING COMPOSITIONAL CHECK

starting_phyla_table <- table(tax_table(Phyloseq_object)[, "phylum"], exclude = NULL)

############### FILTERS

# Filtering obvious taxa that should not be here or that have not been properly classified
filterPhylum <- "unclassified_Bacteria"
filterOrder <- "Chloroplast"
filterFamily <- "Mitochondria"

Phyloseq_filt <- subset_taxa(Phyloseq_object, !is.na(phylum) | phylum != filterPhylum)
Phyloseq_filt <- subset_taxa(Phyloseq_filt, order != filterOrder)
Phyloseq_filt <- subset_taxa(Phyloseq_filt, family != filterFamily)

# Now determine who are the NAs ASVs (from the ASVs fasta file with sequences) and select them
NAs <- subset_taxa(Phyloseq_object, is.na(phylum) | phylum == filterPhylum)
NAs <- taxa_names(NAs)
out <- match(NAs, names(ASVs))
ASVs_NA <- ASVs[out]
rm(NAs, out)

# TODO: maybe this is not totally correct, if you have samples of different types. So maybe it'd be better to separate samples of different types like mosquitoes and water before filtering on prevalence.
# Build prevalence graph for prevalence filtering
prevdf <- apply(X = otu_table(Phyloseq_filt), MARGIN = ifelse(taxa_are_rows(Phyloseq_filt), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(Phyloseq_filt), tax_table(Phyloseq_filt))
prev_graph <- ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(Phyloseq_filt), color = phylum))

# Now selecting the actual good ASVs that have been retained in the dataset for the analysis (no Chloroplast, Mitochondria or NA)
good <- taxa_names(Phyloseq_filt)
keep <- match(good, names(ASVs))
ASVs_good <- ASVs[keep]
rm(good, keep)

# Make VST Phyloseq comparable to this filtered one
Phyloseq_filt_vst <- prune_taxa(taxa_names(Phyloseq_filt), Phyloseq_object_vst)

############### OUTPUTS

setwd(output_path)

png("prevalence_graph.png")

prev_graph + geom_hline(yintercept = 0.02, alpha = 0.2, linetype = 2) + geom_point(size = 1, alpha = 0.7) + scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence / NSamples") + facet_wrap(~phylum) + theme(legend.position = "none")

dev.off()

# Also, write a fasta file containing all the sequences that have not been identified
write.fasta(sequences=ASVs_NA, names=names(ASVs_NA), file.out="ASVs_NA.fasta")
# And a fasta containing the sequences identified by SILVA
write.fasta(sequences=ASVs_good, names=names(ASVs_good), file.out="ASVs_good.fasta")

write.table(starting_phyla_table, file="starting_phyla_table.tsv")

save(Phyloseq_object, Phyloseq_filt, Phyloseq_filt_vst, file="Phyloseq.RData")
save.image(file = "phyl.RData")
