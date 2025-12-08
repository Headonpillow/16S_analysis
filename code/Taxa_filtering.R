############### ENVIRONMENT SETUP (loading libraries and all the snakemake params)

# Libraries
suppressPackageStartupMessages({
    library(phyloseq)
    library(dplyr)
    library(forcats)
    library(ggplot2)
    library(DESeq2)
    library(seqinr)
    library(vegan)
})

main <- function(input_paths = list(), output_paths = list(), params = list()) {
    # Map Snakemake inputs
    samples_file    <- input_paths[["samples"]]    # intermediate/trimmed/samples.txt
    asv_counts_file <- input_paths[["counts"]]     # results/asv/ASVs_counts.tsv
    asv_tax_file    <- input_paths[["tax"]]        # results/asv/ASVs_taxonomy.tsv
    asv_fasta_file  <- input_paths[["fa"]]         # results/asv/ASVs.fa
    metadata_file   <- input_paths[["metadata"]]   # data/meta/metadata.tsv

    # Map Snakemake outputs
    starting_phyla_file <- output_paths[["starting_phyla"]]  # results/phyloseq/starting_phyla_table.tsv
    prevalence_file     <- output_paths[["prevalence"]]      # results/phyloseq/prevalence_graph.png
    phyloseq_file       <- output_paths[["phyloseq"]]        # results/phyloseq/Phyloseq.RData
    asv_good_file       <- output_paths[["asv_good"]]        # results/phyloseq/ASVs_good.fasta

    # Basic checks
    if (any(vapply(list(samples_file, asv_counts_file, asv_tax_file,
                        asv_fasta_file, metadata_file),
                   is.null, logical(1)))) {
        stop("Missing one or more required input_paths entries: expect names 'samples', 'counts', 'tax', 'fa', 'metadata'.")
    }

    if (any(vapply(list(starting_phyla_file, prevalence_file,
                        phyloseq_file, asv_good_file),
                   is.null, logical(1)))) {
        stop("Missing one or more required output_paths entries: expect names 'starting_phyla', 'prevalence', 'phyloseq', 'asv_good'.")
    }
    
    # Resolve paths
    # Use the directory of the phyloseq file as the general output dir (it should have been created from previous scripts)
    output_dir <- dirname(phyloseq_file)
    # Ensure it exists before proceeding
    if (!dir.exists(output_dir)) {
        stop(sprintf("Output directory not found: %s", output_dir))
    }

    ############### PHYLOSEQ OBJECT CREATION
    # Start with loading sample names
    samples <- scan(file = samples_file, what = "character")

    # Now loading the ASVs count_table and the ASVs tax_table
    count_tab <- read.table(asv_counts_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t", fill = TRUE)
    # Replace colnames of count_table with the actual sample names
    colnames(count_tab) <- samples
    tax_tab <- as.matrix(read.table(asv_tax_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t"), mode = "character")
    # Also, read the ASVs sequences for further processing and for selecting the NAs to blast
    ASVs <- read.fasta(asv_fasta_file)

    # Load the metadata and check which ones match with the sample names of the files (keep only those)
    sample_info_tab <- read.table(metadata_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t", fill = TRUE)
    sample_info_tab <- sample_info_tab[rownames(sample_info_tab) %in% samples, , drop = FALSE]
    missing_samples <- setdiff(colnames(count_tab), rownames(sample_info_tab))
    if (length(missing_samples) > 0) {
        stop(paste("Missing metadata for samples:", paste(missing_samples, collapse = ", ")))
    }
    sample_info_tab <- sample_info_tab[colnames(count_tab), , drop = FALSE]
    if (!identical(colnames(count_tab), rownames(sample_info_tab))) {
        stop("Sample names in metadata and count table do not match after filtering.")
    }

    # Finally put all the informations together in a Phyloseq object
    Phyloseq_object <- phyloseq(otu_table(count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab))

    # Filter all the samples with 0 counts
    zero_counts <- colSums(otu_table(Phyloseq_object)) == 0
    Phyloseq_object <- prune_samples(!zero_counts, Phyloseq_object)
    condition <- all(zero_counts == FALSE)

    if (condition == TRUE) {
        cat("\n")
        cat("No samples were removed.")
        cat("\n")
    } else {
        cat("\n")
        filtered <- which(zero_counts == TRUE)
        cat(paste("WARN: sample", names(zero_counts[filtered]), "have been filtered out because it had 0 counts."))
        cat("\n")
    }

    rm(condition)

    ############### DESEQ2
    # Retrieve data from the filtered Phyloseq in data.frame format:
    sample_info_tab <- data.frame(sample_data(Phyloseq_object))
    count_tab <- as.data.frame(otu_table(Phyloseq_object))

    # The formula has been set to "formula ~1" because this equals to no formula, since we are not using this data for differential expression.
    deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~1)
    # If sparse table causes issues, use poscounts
    deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
    deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
    vst_trans_count_tab <- assay(deseq_counts_vst)

    # Now create a phyloseq object with normalized counts
    Phyloseq_object_vst <- phyloseq(otu_table(vst_trans_count_tab, taxa_are_rows = TRUE), sample_data(sample_info_tab), tax_table(tax_tab))

    ############### STARTING COMPOSITIONAL CHECK
    starting_phyla_table <- table(tax_table(Phyloseq_object)[, "phylum"], exclude = NULL)

    ############### FILTERS
    # Filtering obvious taxa that should not be here or that have not been properly classified
    Phyloseq_filt <- subset_taxa(Phyloseq_object, !is.na(phylum) & phylum != "unclassified_Bacteria")
    Phyloseq_filt <- subset_taxa(Phyloseq_filt, order != "Chloroplast")
    Phyloseq_filt <- subset_taxa(Phyloseq_filt, family != "Mitochondria")

    # Now determine who are the NAs ASVs (from the ASVs fasta file with sequences) and select them
    tax_df <- data.frame(tax_table(Phyloseq_object))
    condition <- any(is.na(tax_df$phylum) | tax_df$phylum == "unclassified_Bacteria")

    if (condition == TRUE) {
        cat("\n")
        cat("Some ASVs were not identified, or not bacteria, they have been removed and stored in ASVs_NA")
        cat("\n")
        NAs <- subset_taxa(Phyloseq_object, is.na(phylum) | phylum == "unclassified_Bacteria")
        NAs <- taxa_names(NAs)
        out <- match(NAs, names(ASVs))
        ASVs_NA <- ASVs[out]
        rm(NAs, out)
        # write a fasta file containing all the sequences that have not been identified
        write.fasta(sequences = ASVs_NA, names = names(ASVs_NA), file.out = file.path(output_dir, "ASVs_NA.fasta"))
    }

    # Build prevalence graph for prevalence filtering
    prevdf <- apply(X = otu_table(Phyloseq_filt), MARGIN = ifelse(taxa_are_rows(Phyloseq_filt), yes = 1, no = 2), FUN = function(x) { sum(x > 0) })
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
    png(prevalence_file)
    print(prev_graph + geom_hline(yintercept = 0.02, alpha = 0.2, linetype = 2) + geom_point(size = 1, alpha = 0.7) + scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence / NSamples") + facet_wrap(~phylum) + theme(legend.position = "none"))
    dev.off()

    # And a fasta containing the sequences identified by SILVA
    write.fasta(sequences = ASVs_good, names = names(ASVs_good), file.out = asv_good_file)

    write.table(starting_phyla_table, file = starting_phyla_file, quote = FALSE, sep = "\t", col.names = NA)

    save(Phyloseq_object, Phyloseq_filt, Phyloseq_filt_vst, file = phyloseq_file)
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
