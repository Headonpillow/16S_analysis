############### ENVIRONMENT SETUP (loading libraries and all the snakemake params)

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(ape)
  library(vegan)
})

main <- function(input_paths = list(), output_paths = list(), params = list()) {
  # Map Snakemake inputs
  phyloseq_file <- input_paths[["phyloseq"]]    # results/phyloseq/Phyloseq.RData

  # Map Snakemake params and outputs
  in_dir_param <- params[["in_dir"]]
  out_dir_param <- params[["out_dir"]]

  if (is.null(phyloseq_file)) stop("Missing required input: 'phyloseq'")

  in_path <- if (!is.null(in_dir_param)) in_dir_param else dirname(phyloseq_file)
  out_path <- if (!is.null(out_dir_param)) out_dir_param else dirname(phyloseq_file)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  ############### PHYLOSEQ ANALYSIS
  # Load RData containing Phyloseq objects
  load(phyloseq_file)

  # Also, removing columns in the metadata if one column is all NAs
  info <- sample_data(Phyloseq_filt)
  rm <- names(info)[sapply(info, function(x) sum(is.na(x)) == length(x))]
  info <- info[, !colnames(info) %in% rm]
  sample_data(Phyloseq_filt) <- info
  sample_data(Phyloseq_filt_vst) <- info
  rm(info, rm)

  # Save them again
  save(Phyloseq_filt, Phyloseq_filt_vst, file = "Phyloseq.RData")

  ############### BETA DIVERSITY PRELIMINARY FIGURES
  metadata_fields <- colnames(data.frame(sample_data(Phyloseq_filt_vst)))
  distances <- c("euclidean", "manhattan")

  for (distance in distances) {
    name <- paste(distance, "OD", sep = ".")
    assign(name, ordinate(Phyloseq_filt_vst, "PCoA", distance))
    rm(name, distance)
  }

  ordinations <- list(euclidean.OD, manhattan.OD)

  plot_list <- list()
  i <- 0
  name <- 0

  for (metadata in metadata_fields) {
    for (ordination in ordinations) {
      i <- i + 1
      name <- name + 1
      plot <- plot_ordination(Phyloseq_filt_vst, ordination, type = "samples", color = metadata)
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
    i <- i + 1
    file_name = paste("plot_", i, ".tiff", sep = "")
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
