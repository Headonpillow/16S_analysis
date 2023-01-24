############### ENVIRONMENT SETUP (loading libraries and all the snakemake params)

# Libraries
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(vegan))

# Snakemake paths and params
script_path <- getwd()
IO_dir <- snakemake@params[["IO_dir"]]
IO_path <- paste0(script_path, "/", IO_dir)

############### PHYLOSEQ ANALYSIS

# Start with loading RData containing Phyloseq objects and the tree
setwd(IO_path)
load(file="Phyloseq.RData")
tree <- read.tree(file="ASV_alignment.mafft.treefile")

tree(Phyloseq_filt) <- tree
Phyloseq_filt

# # Now loading the ASVs count_table and the ASVs tax_table as well as the phylogenetic tree
# setwd(asv_files_path)
# count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t")
# # Replace colnames of count_table with the actual sample names
# colnames(count_tab) <- samples
# tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T, row.names=1, check.names=F, sep="\t"))
# tree <- read.tree(file="ASV_tree.nwk")
# # Also, read the ASVs sequences for further processing and for selecting the NAs to blast
# ASVs <- read.fasta("ASVs.fa")
#
# # Load the metadata and check which ones match with the sample names of the files (keep only those)
# setwd(metadata_path)
# sample_info_tab <- read.table("metadata.tsv", header=T, row.names=1, check.names=F, sep="\t")
# sample_info_tab <- filter(sample_info_tab, rownames(sample_info_tab) %in% samples)
# order <- match(rownames(sample_info_tab), colnames(count_tab))
# sample_info_tab <- arrange(sample_info_tab, order)
#
# # Finally put all the informations together in a Phyloseq object
# Phyloseq_object <- phyloseq(otu_table(count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab), tree)
#
# ############### DESEQ2
#
# deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~breeding_site)
# # Since getting an error with the dataset mosquitoes+water:
# # "Error in estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# # gene contains at least one zero, cannot compute log geometric means"
# # It's because the table is quite sparse with many zeroes, so I needed to add the
# # Next command to deal with it.
# # If there is not such problem, next line can be commented.
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# vst_trans_count_tab <- assay(deseq_counts_vst)
#
# # Now create a phyloseq object with normalized counts
# Phyloseq_object_vst <- phyloseq(otu_table(vst_trans_count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab), tree)
#
# ############### STARTING COMPOSITIONAL CHECK
#
# starting_phyla_table <- table(tax_table(Phyloseq_object)[, "phylum"], exclude = NULL)
#
# ############### FILTERS
#
# # Filtering obvious taxa that should not be here or that have not been properly classified
# filterPhylum <- "unclassified_Bacteria"
# filterOrder <- "Chloroplast"
# filterFamily <- "Mitochondria"
#
# Phyloseq_filt <- subset_taxa(Phyloseq_object, !is.na(phylum) | phylum != filterPhylum)
# Phyloseq_filt <- subset_taxa(Phyloseq_filt, order != filterOrder)
# Phyloseq_filt <- subset_taxa(Phyloseq_filt, family != filterFamily)
#
# # Now determine who are the NAs ASVs (from the ASVs fasta file with sequences) and select them
# NAs <- subset_taxa(Phyloseq_object, is.na(phylum) | phylum == filterPhylum)
# NAs <- taxa_names(NAs)
# keep <- match(NAs, names(ASVs))
# ASVs_NA <- ASVs[keep]
#
# # Build prevalence graph for prevalence filtering
# prevdf <- apply(X = otu_table(Phyloseq_filt), MARGIN = ifelse(taxa_are_rows(Phyloseq_filt), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(Phyloseq_filt), tax_table(Phyloseq_filt))
# prev_graph <- ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(Phyloseq_filt), color = phylum))
#
# # Filter by prevalence setting a threshold of 3%
# # Basically this kind of filtering is supposed to improve the downstream pipeline removing contaminants.
# # The filter removes taxa that are present in less than 3% of the samples. (so very rare sequences)
# prevalenceThreshold <- 0.03 * nsamples(Phyloseq_filt)
# keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
# Phyloseq_filt <- prune_taxa(keepTaxa, Phyloseq_filt)
#
# # Make VST Phyloseq comparable to this filtered one
# Phyloseq_filt_vst <- prune_taxa(taxa_names(Phyloseq_filt), Phyloseq_object_vst)
#
# ############### BETA DIVERSITY PRELIMINARY FIGURES
#
#
# metadata_fields <- as.list(sample_info_tab[0,])
# distances <- c("euclidean", "manhattan", "unifrac", "wunifrac")
#
# calculate_dist_matrix <- function(physeq, dist_measures){
#   a <- list()
#   for (distance in dist_measures){
#     name <- paste(distance, "OD", sep=".")
#     append(a, assign(name, ordinate(physeq, "PCoA", distance)))
#   }
#   print(matrix_list)
# }
#
# calculate_dist_matrix(Phyloseq_filt_vst, distances)
#
# print(matrix_list)
#
#
# #### But this would work only for one specific thing selecting individuals etc.
# #### I need a way to generalize this and get combination of useful metadata, could
# #### As well read the whole metadata row or levels and do combinations for all of them??
#
# # # Select individuals only
# # individual_vst <- prune_samples(Phyloseq_filt_vst, sample_type=="Individual")
# #
# # adonis_df_ind_vst <- as.data.frame(t(otu_table(individual_vst)))
# # adonis_info_ind_vst <- sample_data(individual_vst)
# # individual_vst.dist <- vegdist(adonis_df_ind_vst, method="euclidean")
# #
# # # Making the info df as data.frame
# # adonis_info_ind_vst <- as.matrix(sample_data(individual_vst))
# # adonis_info_ind_vst <- as.data.frame(adonis_info_ind_vst)
# #
# # individual_vst.div <- adonis2(adonis_df_ind_vst ~breeding_site*breeding_site_type, data=adonis_info_ind_vst, permutations=999, method="euclidean")
# # dispersion <- betadisper(individual_vst.dist, group=adonis_info_ind_vst$breeding_site)
# # permutest(dispersion)
# # plot(dispersion, hull=FALSE, ellipse=TRUE)
#
# ############### OUTPUTS
#
# setwd(output_path)
#
# png("prevalence_graph.png")
#
# prev_graph + geom_hline(yintercept = 0.02, alpha = 0.2, linetype = 2) + geom_point(size = 1, alpha = 0.7) + scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence / NSamples") + facet_wrap(~phylum) + theme(legend.position = "none")
#
# dev.off()
#
# # Also, write a fasta file containing all the sequences that have not been identified
# write.fasta(sequences=ASVs_NA, names=names(ASVs_NA), file.out="ASVs_NA.fasta")
#
# write.table(starting_phyla_table, file="starting_phyla_table.tsv")
#
# save(Phyloseq_filt, Phyloseq_filt_vst, file="Phyloseq.RData")
