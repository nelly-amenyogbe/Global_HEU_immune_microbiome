###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Microbiome richness and diversity
###################

# In this script, we will compute the Shannon Diversity for all individuals and perform the Wilcoxon Test to determine whether either richness or diversity differ between HEU and HUU children in each study site separately.

# load packages
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)

# load data
sample.cols <- read.csv("HEU_manuscript_analysis/Rdata/luminex/heu_pca_colours.csv")
physeq <- readRDS("HEU_manuscript_analysis/Rdata/microbiome/gc_heu_phyloseq.rds")

# Subset data by site ####
cad <- subset_samples(physeq, Site == "Canada")
cad <- prune_taxa(taxa_sums(cad) > 0, cad)

blg <- subset_samples(physeq, Site == "Belgium")
blg <- prune_taxa(taxa_sums(blg) > 0, blg)

saf <- subset_samples(physeq, Site == "South Africa")
saf <- prune_taxa(taxa_sums(saf) > 0, saf)

sample_data(physeq)$site.exposure <- paste(sample_data(physeq)$Site, sample_data(physeq)$Exposure)

# Determine sample sequencing depth ####

# prepare data
# get sample depth
meta <- data.frame(sample_data(physeq))
depth <- colSums(otu_table(physeq))
meta$depth <- depth

p.depth <- ggplot(meta, aes(x = Site, y = depth, color = Exposure)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right") +
  labs(x = "", y = "Number of reads")

p.depth

# conclude: There is one sample within the Canadian cohort with much lower sequencing depth than the remainder of samples.  We will remove this sample from the dataset, and for the remainder of the samples, rarefy the data prior to calulating diversity estimates.

# Rarefy data ####
gc.rar <- prune_samples(sample_sums(physeq) > 20000, physeq)
gc.rar <- rarefy_even_depth(gc.rar, rngseed = 711)

rar.sub <- colSums(otu_table(gc.rar))
unique(rar.sub) # rarified to 33113 reads per sample

# Calculate Shannon Diversity ####

shan.rich <- plot_richness(gc.rar, measures = "Shannon")
shan.rich.data <- shan.rich$data

# Plot diversity by cohort ####
# Change control to HUU

# prepare graphing aesthetics
sample.cols$group <- gsub("Control", "HUU", sample.cols$group)
sample.cols$site.heu <- paste(sample.cols$site, sample.cols$group)
p.cols <- as.character(sample.cols$colour)
names(p.cols) <- sample.cols$site.heu

# prepare factor levels in richness data frame
shan.rich.data$Exposure <- gsub("Control", "HUU", shan.rich.data$Exposure)
shan.rich.data$site.exposure <- paste(shan.rich.data$Site, shan.rich.data$Exposure)
shan.rich.data$site.exposure <- factor(shan.rich.data$site.exposure, levels = as.character(sample.cols$site.heu))

# plot
p.shannon <- ggplot(shan.rich.data, aes(x = Site, y = value, fill = site.exposure)) +
  geom_boxplot(outlier.size = NA, color = "black") +
  geom_point(shape = 21, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  scale_fill_manual("", values = p.cols) + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.ticks.x = element_blank(),
        legend.position = "right") +
  labs(x = "", y = "Shannon Index")

p.shannon

#ggsave("HEU_manuscript_analysis/figures/microbiome/heu_shannon_boxplots.pdf", device = "pdf", dpi = 300, width = 4, height = 3)

#### Shannon Diversity Statistics ####
sites <- as.character(unique(shan.rich.data$Site))

rich.wilcox <- ldply(as.list(sites), function(i){
  df <- filter(shan.rich.data, Site == i)
  heu <- filter(df, Exposure == "HEU")
  ue <- filter(df, Exposure == "HUU")
  test <- wilcox.test(heu$value, ue$value)
  test$p.value
})

rich.wilcox$Site <- sites
rich.wilcox

# Conclude:  No significant differences in Shannon Diversity between HEU and HUU children within any site
 
# END ####