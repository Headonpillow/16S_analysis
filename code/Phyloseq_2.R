# PREVALENCE FILTERING

# Here you first extract an OTU table from the phyloseq object Tanzania (for getting read counts for sample, for OTU), then ceck where are the taxa in the phyloseq object Tanzania, they're on columns.
# Then, you apply on the columns a function to sum all the times the OTU read count was above zero. (counting on samples, so MAX is 74 times in this case)
prevdf <- apply(X = otu_table(Tanzania_filt), MARGIN = ifelse(taxa_are_rows(Tanzania_filt), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# De facto, the resulting object is indeed a vector of integers, one for each OTU. This number are prevalence numbers, or the samples each OTU had a read count above 0 in our dataset.

# Add taxonomy columns (one for each rank) and use the taxa_sums to return the total read count (Abundance) for each OTU (across the entire dataset), then make it a data frame
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(Tanzania_filt), tax_table(Tanzania_filt))

# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.02 * nsamples(Tanzania_filt)
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]

Tanzania_filt <- prune_taxa(keepTaxa, Tanzania_filt)


## ----------------------------------------------------------------------
# Filter Wolbachia
Tanzania_filt <- subset_taxa(Tanzania_filt, genus != "Wolbachia")

# Remove useless samples
samples_to_remove <- c("Tanz_48", "Tanz_30", "Tanz_26", "Tanz_45", "Tanz_51", "Tanz_57")
Tanzania_filt <- subset_samples(Tanzania_filt, !sample_names(Tanzania_filt) %in% samples_to_remove)

# Create two different datasets for swabs and mosquitoes
Tanzania_swabs <- subset_samples(Tanzania_filt, category=="swabs")
Tanzania_mosquitoes <- subset_samples(Tanzania_filt, category=="mosquitoes")

# Filter mosquitoes with no taxonomical information from the COI
Tanzania_mosquitoes <- subset_samples(Tanzania_mosquitoes, mos_gen != "")

# change numbers to roman

sample_data <- sample_data(Tanzania_mosquitoes)
sample_data <- data.frame(sample_data)
sample_data <- mutate(sample_data, house=case_when(house == "I" ~ 1,  house == "II" ~ 2, house == "III" ~ 3, house == "IV" ~ 4, house == "V" ~ 5, house == "VI" ~ 6, house == "VII" ~ 7, house == "IX" ~ 9, house == "X" ~ 10))
sample_data$house <- as.factor(sample_data$house)
sample_data(Tanzania_mosquitoes) <- sample_data

Tanzania_mosquitoes_rel <- transform_sample_counts(Tanzania_mosquitoes, function(x) 1E6 * x/sum(x))
Tanzania_mosquitoes_glom <- tax_glom(Tanzania_mosquitoes, taxrank = "genus")
Tanzania_mosquitoes_glom_rel <- transform_sample_counts(Tanzania_mosquitoes_glom, function(x) x/sum(x))
# With this filter you go down from 82 to 71 samples

# Separate them also according to the different genus
Tanzania_culex <- subset_samples(Tanzania_mosquitoes, mos_gen != "Culex")
Tanzania_aedes <- subset_samples(Tanzania_mosquitoes, mos_gen != "Aedes")

# # Also, after relative transformation some of them are NaN
# head(otu_table(Tanzania_mosquitoes_rel)[3])
# head(otu_table(Tanzania_mosquitoes_glom_rel)[3])
#
# # In both cases are 48, 30, 26, 45, 51, 57
#
# enough <- filter(data.frame(otu_table(Tanzania_mosquitoes_glom)[,10]), Tanz_48 > 1)
# sum(enough)
#
# # these samples have in fact abundance of 0 for all the taxa, so they can as well get filtered at an earlier stage


## ----------------------------------------------------------------------
Tanzania_mosquitoes_glom


## ----------------------------------------------------------------------
reads_per_sample <- sample_sums(Tanzania_mosquitoes_glom)
reads_per_sample

sample_data <- data.frame(sample_data(Tanzania_mosquitoes_glom))
sample_data <- mutate(sample_data, reads = reads_per_sample)
sample_data$reads <- as.integer(sample_data$reads)
aggregate_df <- sample_data[c("house", "reads")]
reads_per_house <- aggregate(aggregate_df[-1], by = list(aggregate_df$house), FUN = sum, na.rm = T)
reads_per_house

sample_per_house <- table(sample_data$house)
sample_per_house


## ----------------------------------------------------------------------
samples_to_keep <- c("Tanz_62", "Tanz_78", "Tanz_81", "Tanz_86", "Tanz_90", "Tanz_82", "Tanz_80", "Tanz_92", "Tanz_63", "Tanz_68")
Tanzania_mosquitoes <- subset_samples(Tanzania_mosquitoes, sample_names(Tanzania_mosquitoes) %in% samples_to_keep)
Tanzania_mosquitoes_glom <- subset_samples(Tanzania_mosquitoes_glom, sample_names(Tanzania_mosquitoes_glom) %in% samples_to_keep)
Tanzania_mosquitoes_rel <- subset_samples(Tanzania_mosquitoes_rel, sample_names(Tanzania_mosquitoes_rel) %in% samples_to_keep)


## ----------------------------------------------------------------------
# sample_sums(Tanzania_mosquitoes)
#
# a <- prune_samples(sample_sums(Tanzania_mosquitoes) > 10, Tanzania_mosquitoes)


## ----------------------------------------------------------------------
# diagdds = phyloseq_to_deseq2(a, ~ house)
#
# # calculate geometric means prior to estimate size factors
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# geoMeans = apply(counts(diagdds), 1, gm_mean)
# diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans, type = "iterate")
# diagdds = DESeq(diagdds, fitType="local")
#
# res = results(diagdds)
# res = res[order(res$padj, na.last=NA), ]
# alpha = 0.01
# sigtab = res[(res$padj < alpha), ]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(a)[rownames(sigtab), ], "matrix"))
# head(sigtab)



## ----------------------------------------------------------------------
Tanzania_nine_swabs <- merge_phyloseq(Tanzania_mosquitoes, Tanzania_swabs)

Tanzania_nine_swabs <- rarefy_even_depth(Tanzania_nine_swabs, sample.size = 1000, rngseed = 10)
Tanzania_mosquitoes <- rarefy_even_depth(Tanzania_mosquitoes, sample.size = 1000, rngseed = 10)


## ----------------------------------------------------------------------
sample_data <- sample_data(Tanzania_mosquitoes)
sample_data <- data.frame(sample_data)
sample_data <- mutate(sample_data, spec=case_when(mos_spec == "egypti" ~ "Ae. aegypti",  mos_spec == "hybrid" ~ "Ae. hybrid", mos_spec == "pipiens" ~ "Cx. pipiens", mos_spec == "quinquefasciatus" ~ "Cx. quinquefasciatus", TRUE ~ as.character(mos_spec)))

sample_data(Tanzania_mosquitoes) <- sample_data

#Tanzania_aedes_rel <-  subset_samples(Tanzania_mosquitoes_rel, mos_gen != "Culex")
#Tanzania_culex_rel <- subset_samples(Tanzania_mosquitoes_rel, mos_gen != "Aedes")

Tanzania_mosquitoes.ord <- ordinate(Tanzania_mosquitoes, "PCoA", "bray")
#Tanzania_aedes.ord <- ordinate(Tanzania_aedes_rel, "PCoA", "bray")
#Tanzania_culex.ord <- ordinate(Tanzania_culex_rel, "PCoA", "bray")

p1 <- plot_ordination(Tanzania_mosquitoes, Tanzania_mosquitoes.ord, type="samples", color="house") + geom_point(size=3) + labs(title = "PCoA of Houses (Bray)", color = "House") + theme(plot.title = element_text(hjust = 0.5))
#p2 <- plot_ordination(Tanzania_aedes, Tanzania_aedes.ord, type="samples", color="house") + geom_point(size=3) + labs(title = "PCoA of houses - Aedes", color = "House") + theme(plot.title = element_text(hjust = 0.5))
#p3 <- plot_ordination(Tanzania_culex, Tanzania_culex.ord, type="samples", color="house") + geom_point(size=3) + labs(title = "PCoA of houses - Culex", color = "House") + theme(plot.title = element_text(hjust = 0.5))

p1
#p2
#p3

## ----------------------------------------------------------------------
# Tanzania_mosquitoes.ord <- ordinate(Tanzania_mosquitoes_rel, "NMDS", "bray", weighted = TRUE)
# Tanzania_aedes.ord <- ordinate(Tanzania_aedes_rel, "NMDS", "bray")
# Tanzania_culex.ord <- ordinate(Tanzania_culex_rel, "NMDS", "bray")
#
# p4 = plot_ordination(Tanzania_mosquitoes_rel, Tanzania_mosquitoes.ord, type="samples", color="spec") + geom_point(size=3) + labs(title = "NMDS of species", color = "Specie") + theme(plot.title = element_text(hjust = 0.5))
# p5 = plot_ordination(Tanzania_aedes_rel, Tanzania_aedes.ord, type="samples", color="spec") + geom_point(size=3) + labs(title = "NMDS of species - Aedes", color = "Specie") + theme(plot.title = element_text(hjust = 0.5))
# p6= plot_ordination(Tanzania_culex_rel, Tanzania_culex.ord, type="samples", color="spec") + geom_point(size=3) + labs(title = "NMDS of species - Culex", color = "Specie") + theme(plot.title = element_text(hjust = 0.5))
#
# p4
# p5
# p6


## ----------------------------------------------------------------------
Tanzania_mosquitoes.ord <- ordinate(Tanzania_mosquitoes, "PCoA", "unifrac", weighted = TRUE)
#Tanzania_aedes.ord <- ordinate(Tanzania_aedes_rel, "PCoA", "unifrac")
#Tanzania_culex.ord <- ordinate(Tanzania_culex_rel, "PCoA", "unifrac")

p7 <- plot_ordination(Tanzania_mosquitoes, Tanzania_mosquitoes.ord, type="samples", color="house") + geom_point(size=3) + labs(title = "PCoA of Houses (Unifrac)", color = "House") + theme(plot.title = element_text(hjust = 0.5))
#p8 <- plot_ordination(Tanzania_aedes_rel, Tanzania_aedes.ord, type="samples", color="spec") + geom_point(size=3) + labs(title = "PCoA of species - Aedes", color = "Specie") + theme(plot.title = element_text(hjust = 0.5))
#p9 <- plot_ordination(Tanzania_culex_rel, Tanzania_culex.ord, type="samples", color="spec") + geom_point(size=3) + labs(title = "PCoA of species - Culex", color = "Specie") + theme(plot.title = element_text(hjust = 0.5))

p7
#p8
#p9


## ----------------------------------------------------------------------
#Indicator species

ASV_df <- data.frame(t(otu_table(Tanzania_mosquitoes)))
meta_df <- data.frame(sample_data(Tanzania_mosquitoes))
indicator_df <- cbind(house=meta_df$house, ASV_df)
indicator_df <- cbind(mos_gen=meta_df$mos_gen, indicator_df)
indicator_df <- cbind(sample=rownames(indicator_df), indicator_df)

abund <- indicator_df[,4:ncol(indicator_df)]
house <- indicator_df$house
inv <- multipatt(abund, house, func = "r.g", control = how(nperm=1000))


## ----------------------------------------------------------------------
#### temporary chunk for compositional profiles (Phylum is the good one)

Culex_phylum <- tax_glom(Tanzania_culex, taxrank = "phylum")
Culex_phylum_rel <- transform_sample_counts(Culex_phylum, function(x) x / sum(x))
Culex_class <- tax_glom(Tanzania_culex, taxrank = "class")
Culex_class_rel <- transform_sample_counts(Culex_class, function(x) x / sum(x))
Culex_order <- tax_glom(Tanzania_culex, taxrank = "order")
Culex_order_rel <- transform_sample_counts(Culex_order, function(x) x / sum(x))

Culex_phylum_long <- psmelt(Culex_phylum_rel)
Culex_class_long <- psmelt(Culex_class_rel)
Culex_order_long <- psmelt(Culex_order_rel)

Culex_phylum_long <- mutate(Culex_phylum_long, phylum=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(phylum)))
Culex_class_long <- mutate(Culex_class_long, phylum=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(phylum)), class=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(class)))
Culex_order_long <- mutate(Culex_order_long, phylum=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(phylum)), class=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(class)), order=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(order)))

# Here is to set colors and the right factor orders for phyla
Culex_phyla_col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#000000")

nb.class <- 14
mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.class)
Culex_class_col <- append(mycolors, "#000000")

Culex_phylum_long$phylum <- as.factor(Culex_phylum_long$phylum)
Culex_phylum_long$phylum <- factor(Culex_phylum_long$phylum, levels = c("Acidobacteriota", "Actinobacteriota","Bacteroidota", "Campilobacterota", "Firmicutes", "Myxococcota", "Patescibacteria", "Proteobacteria", "Other"))
Culex_class_long$class <- as.factor(Culex_class_long$class)
Culex_class_long$class <- factor(Culex_class_long$class, levels = c("Acidimicrobiia", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Blastocatellia", "Campylobacteria", "Clostridia", "Gammaproteobacteria", "Myxococcia", "Negativicutes", "Rubrobacteria", "Saccharimonadia", "Thermoleophilia", "Other"))

Culex_phylum_plot <- ggplot(Culex_phylum_long, aes(x = house, y = Abundance, fill = phylum, color = phylum)) +
  geom_col(position = "fill") +
  theme(axis.text.x=element_text(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance") +
  scale_color_manual(values = Culex_phyla_col) +
  scale_fill_manual(values = Culex_phyla_col) +
  labs(title = "Culex phyla") +
  theme(plot.title = element_text(hjust = 0.5))

Culex_class_plot <- ggplot(Culex_class_long, aes(x = house, y = Abundance, fill = class, color = class)) +
  geom_col(position = "fill") +
  theme(axis.text.x=element_text(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance") +
  scale_color_manual(values = Culex_class_col) +
  scale_fill_manual(values = Culex_class_col) +
  labs(title = "Culex class") +
  theme(plot.title = element_text(hjust = 0.5))

Culex_order_plot <- ggplot(Culex_order_long, aes(x = house, y = Abundance, fill = order)) +
  geom_col(position = "fill") +
  theme(axis.text.x=element_text(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance")

Aedes_phylum <- tax_glom(Tanzania_aedes, taxrank = "phylum")
Aedes_phylum_rel <- transform_sample_counts(Aedes_phylum, function(x) x / sum(x))
Aedes_class <- tax_glom(Tanzania_aedes, taxrank = "class")
Aedes_class_rel <- transform_sample_counts(Aedes_class, function(x) x / sum(x))
Aedes_order <- tax_glom(Tanzania_aedes, taxrank = "order")
Aedes_order_rel <- transform_sample_counts(Aedes_order, function(x) x / sum(x))

Aedes_phylum_long <- psmelt(Aedes_phylum_rel)
Aedes_class_long <- psmelt(Aedes_class_rel)
Aedes_order_long <- psmelt(Aedes_order_rel)

Aedes_phylum_long <- filter(Aedes_phylum_long, house != 8)
Aedes_phylum_long <- filter(Aedes_phylum_long, ! is.na(phylum))

Aedes_phylum_long <- mutate(Aedes_phylum_long, phylum=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(phylum)))
Aedes_class_long <- mutate(Aedes_class_long, phylum=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(phylum)), class=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(class)))
Aedes_order_long <- mutate(Aedes_order_long, phylum=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(phylum)), class=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(class)), order=case_when(Abundance < 0.001 | is.na(Abundance) ~ "Other", TRUE ~ as.character(order)))

# Here is to set colors and the right factor orders for phyla
Aedes_phyla_col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99", "#33a02c", "#e31a1c", "#fdbf6f", "#ff7f00", "#000000")
Aedes_phylum_long$phylum <- as.factor(Aedes_phylum_long$phylum)
Aedes_phylum_long$phylum <- factor(Aedes_phylum_long$phylum, levels = c("Acidobacteriota", "Actinobacteriota","Bacteroidota", "Firmicutes", "Fusobacteriota", "Myxococcota", "Patescibacteria", "Plantomycetota", "Proteobacteria", "Other"))
Aedes_class_long$class <- as.factor(Aedes_class_long$class)
Aedes_class_long$class <- factor(Aedes_class_long$class, levels = c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Clostridia", "Coriobacteriia", "Fusobacteriia", "Gammaproteobacteria", "Myxococcia", "Negativicutes", "Planctomycetes", "Polyangia", "Rubrobacteria", "Saccharimonadia", "Thermoleophilia", "Other"))
Aedes_class_col <- c("#7FC97F", "#A0BAAC", "#C2AFCE", "#E4B9A3", "#FDC988", "#658DAA", "#FEEB93", "#D1DD9E", "#704BA0", "#D31286", "#DD2456", "#D23740", "#B95B1C", "#C2541E", "#95603B", "#666666", "#000000")

Aedes_phylum_plot <- ggplot(Aedes_phylum_long, aes(x = house, y = Abundance, fill = phylum, color = phylum)) +
  geom_col(position = "fill") +
  theme(axis.text.x=element_text(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance") +
  scale_color_manual(values = Aedes_phyla_col) +
  scale_fill_manual(values = Aedes_phyla_col) +
  labs(title = "Aedes phyla") +
  theme(plot.title = element_text(hjust = 0.5))

Aedes_class_long <- filter(Aedes_class_long, house != 8)
Aedes_class_plot <- ggplot(Aedes_class_long, aes(x = house, y = Abundance, fill = class, color=class)) +
  geom_col(position = "fill") +
  theme(axis.text.x=element_text(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance") +
  scale_color_manual(values = Aedes_class_col) +
  scale_fill_manual(values = Aedes_class_col) +
  labs(title = "Aedes class") +
  theme(plot.title = element_text(hjust = 0.5))

Aedes_order_plot <- ggplot(Aedes_order_long, aes(x = house, y = Abundance, fill = order)) +
  geom_col(position = "fill") +
  theme(axis.text.x=element_text(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y="Relative abundance")


## ----------------------------------------------------------------------
# plot the mosquitoes according to the coordinates to see if there are any groups

ggplot(data=Tanzania_mosquitoes_rel@sam_data, mapping = aes_string(x = "longitude", y="latitude", color="latitude")) + scale_color_viridis(option="magma") + geom_point(size = 1, position = position_jitter(w=0.0008, h=0.0008)) + geom_violin(fill=NA)


## ----------------------------------------------------------------------
# creating a new variable for North, Center and South

sample_data(Tanzania_mosquitoes_rel)$latitude_binned <- cut(sample_data(Tanzania_mosquitoes_rel)$latitude, breaks = c(-7.8, -7.79, -7.783, -7.775))
levels(sample_data(Tanzania_mosquitoes_rel)$latitude_binned) <- list(South="(-7.8,-7.79]", Center="(-7.79,-7.783]", North="(-7.783,-7.775]")