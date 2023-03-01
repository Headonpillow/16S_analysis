############### ENVIRONMENT SETUP (loading libraries and all the snakemake params)

# Libraries
# suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(phyloseq))
# suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(vegan))

# Snakemake paths and params
script_path <- getwd()
in_dir <- snakemake@params[["in_dir"]]
out_dir <- snakemake@params[["out_dir"]]
in_path <- paste0(script_path, "/", in_dir)
out_path <- paste0(script_path, "/", out_dir)

############### PHYLOSEQ ANALYSIS

# Start with loading RData containing Phyloseq objects and the tree
setwd(in_path)
load(file="Phyloseq.RData")
tree <- read_tree(treefile="ASV_alignment.mafft.treefile")

# Adding the tree into the phyloseq objects for the use of phylogenetic distance measures
Phyloseq_filt <- merge_phyloseq(Phyloseq_filt, tree)
Phyloseq_filt_vst <- merge_phyloseq(Phyloseq_filt_vst, tree)

# Save them again
save(Phyloseq_filt, Phyloseq_filt_vst, file="Phyloseq.RData")

############### BETA DIVERSITY PRELIMINARY FIGURES

metadata_fields <- colnames(data.frame(sample_data(Phyloseq_filt_vst)))
distances <- c("euclidean", "manhattan", "unifrac", "wunifrac")

for (distance in distances) {
  name <- paste(distance, "OD", sep=".")
  assign(name, ordinate(Phyloseq_filt_vst, "PCoA", distance))
  rm(name,distance)
}

# TODO: make some of the fields of the metadata factors, because they are still seen as numeric like "date"
# TODO: plots fail in the case that metadata file is just one column, probably a bug in Phyloseq. Column should be at least two, even if one is just samples. Maybe I can just print a warning.

ordinations <- list(euclidean.OD, manhattan.OD, unifrac.OD, wunifrac.OD)

plot_list <- list()
i <- 0
name <- 0

for (metadata in metadata_fields) {
  for (ordination in ordinations) {
    i <- i+1
    name <- name+1
    plot <- plot_ordination(Phyloseq_filt_vst, ordination, type="samples", color = metadata)
    plot <- plot + ggtitle(distances[name])
    plot_list[[i]] <- plot
  }
  name <- 0
}

rm(i, name)

# Change the output path
setwd(out_path)

i <- 0
for (plot in plot_list) {
  i <- i+1
  file_name = paste("plot_", i, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}

# TODO: fix adonis and think about how to divide the output
############### PERMANOVA

# adonis_df_vst <- as.data.frame(t(otu_table(Phyloseq_filt_vst)))
# adonis_info_vst <- sample_data(Phyloseq_filt_vst)
#
# for (distance in distances) {
#   name <- paste(distance, "dist", sep=".")
#   assign(name, vegdist(adonis_df_vst, method=distance))
#   rm(name,distance)
# }
#
# ordinations <- list(euclidean.OD, manhattan.OD, unifrac.OD, wunifrac.OD)
#
# adonis_info_ind_vst <- as.matrix(sample_data(Phyloseq_filt_vst))
# adonis_info_ind_vst <- as.data.frame(adonis_info_vst)

# individual_vst.div <- adonis2(adonis_df_ind_vst ~breeding_site*breeding_site_type, data=adonis_info_ind_vst, permutations=999, method="euclidean")
# dispersion <- betadisper(individual_vst.dist, group=adonis_info_ind_vst$breeding_site)
# permutest(dispersion)
# plot(dispersion, hull=FALSE, ellipse=TRUE)
