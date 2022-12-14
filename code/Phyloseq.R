############### ENVIRONMENT SETUP (loading libraries and all the snakemake params)

# Libraries
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ape))

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

# Now loading the ASVs count_table and the ASVs tax_table as well as the phylogenetic tree
setwd(asv_files_path)
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t")
# Replace colnames of count_table with the actual sample names
colnames(count_tab) <- samples
tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T, row.names=1, check.names=F, sep="\t"))
tree <- read.tree(file="ASV_tree.nwk")

# Load the metadata and check which ones match with the sample names of the files (keep only those)
setwd(metadata_path)
sample_info_tab <- read.table("metadata.csv", header=T, row.names=1, check.names=F, sep=",")
sample_info_tab <- filter(sample_info_tab, rownames(sample_info_tab) %in% samples)
order <- match(rownames(sample_info_tab), colnames(count_tab))
sample_info_tab <- arrange(sample_info_tab, order)

# Finally put all the informations together in a Phyloseq object
Phyloseq_object <- phyloseq(otu_table(count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab), tree)

############### STARTING COMPOSITIONAL CHECK

starting_phyla_table <- table(tax_table(Phyloseq_object)[, "phylum"], exclude = NULL)

############### FILTERS

# Filtering obvious taxa that should not be here or that have not been properly classified
filterPhylum <- "unclassified_Bacteria"
filterOrder <- "Chloroplast"
filterFamily <- "Mitochondria"

Phyloseq_filt <- subset_taxa(Phyloseq_object, !is.na(phylum) & phylum != filterPhylum)
Phyloseq_filt <- subset_taxa(Phyloseq_filt, order != filterOrder)
Phyloseq_filt <- subset_taxa(Phyloseq_filt, family != filterFamily)

# Build prevalence graph for prevalence filtering
prevdf <- apply(X = otu_table(Phyloseq_filt), MARGIN = ifelse(taxa_are_rows(Phyloseq_filt), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(Phyloseq_filt), tax_table(Phyloseq_filt))
prev_graph <- ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(Phyloseq_filt), color = phylum))

# Filter by prevalence setting a threshold of 3%
# Basically this kind of filtering is supposed to improve the downstream pipeline removing contaminants.
# The filter removes taxa that are present in less than 3% of the samples. (so very rare sequences)
prevalenceThreshold <- 0.03 * nsamples(Phyloseq_filt)
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
Phyloseq_filt <- prune_taxa(keepTaxa, Phyloseq_filt)

############### OUTPUTS

setwd(output_path)

png("prevalence_graph.png")

prev_graph + geom_hline(yintercept = 0.02, alpha = 0.2, linetype = 2) + geom_point(size = 1, alpha = 0.7) + scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence / NSamples") + facet_wrap(~phylum) + theme(legend.position = "none")

dev.off()

write.table(starting_phyla_table, file="starting_phyla_table.tsv")

save(Phyloseq_object, file="Phyloseq.RData")