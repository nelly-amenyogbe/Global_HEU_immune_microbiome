###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Microbiome data ordination
###################

# In this script, we will prepare ordinations (Non-metric multidimensional scaling; NMDS) for HEU and HUU children together, and for each study site separately.  We will then perform statistical analyses (PERMANOVA test) to determine the variance explained by HIV exposure on the stool microbiome composition in each study site separately.

# Load packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(plyr)

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

# Ordination of samples ####
# get.nmds.data
# Input: ps = phyloseq object
# output: data frame friendly for plotting via ggplot2

get.nmds.data <- function(ps){
  ord <- ordinate(ps, method = "NMDS", distance = "bray")
  p <- plot_ordination(physeq, ord)
  dat <- p$data
  dat$Exposure <- gsub("Control", "HUU", dat$Exposure)
  dat$Exposure <- factor(dat$Exposure, levels = c("HUU", "HEU"))
  dat$site.exposure <- paste(dat$Site, dat$Exposure)
  return(dat)
}

# All samples together ####

heu.nmds.dat <- get.nmds.data(physeq)

# set aesthetics
# colour
sample.cols$site.heu <- gsub("Control", "HUU", sample.cols$site.heu)

cols <- as.character(sample.cols$colour)
names(cols) <- sample.cols$site.heu

# shape
shapes <- sample.cols$shape
names(shapes) <- sample.cols$site.heu

shapes.f <- c("Canada HEU" = 24, "Canada HUU" = 24, "Belgium HEU" = 21, "Belgium HUU" = 21, "South Africa HEU" = 22, "South Africa HUU" = 22) # this is for shapes where the fill is specified rather than the colour

# set factor levels
heu.nmds.dat$site.exposure <- factor(heu.nmds.dat$site.exposure, levels = c("Belgium HUU", "Belgium HEU", "Canada HUU", "Canada HEU", "South Africa HUU", "South Africa HEU"))

p.heu <- ggplot(heu.nmds.dat, aes(x = NMDS1, y =NMDS2, fill = site.exposure, shape = site.exposure)) +
  geom_point(alpha = 0.9, size = 4, color = "black") +
  theme_classic() +
  scale_shape_manual("", values = shapes.f) +
  scale_fill_manual("", values = cols) + 
  theme(axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom")

p.heu

#ggsave("HEU_manuscript_analysis/figures/microbiome/heu_nmds.pdf", device = "pdf", dpi = 300, width = 6.5, height = 6)

# Ordination by Site ####

# CAD ord ####
cad.dat <- get.nmds.data(cad)

cad.p <- ggplot(cad.dat, aes(x = NMDS1, y =NMDS2, fill = site.exposure, shape = site.exposure)) +
  geom_point(alpha = 0.9, size = 4, shape = 24) +
  theme_classic() +
  scale_shape_manual("", values = shapes[3:4]) +
  scale_fill_manual("", values = cols[3:4]) + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.position = "none")

cad.p

#ggsave("HEU_manuscript_analysis/figures/microbiome/cad_heu_nmds.pdf", device = "pdf", dpi = 300, width = 5, height = 4)

# BLG ord ####
blg.dat <- get.nmds.data(blg)

blg.p <- ggplot(blg.dat, aes(x = NMDS1, y =NMDS2, fill = site.exposure, shape = site.exposure)) +
  geom_point(alpha = 0.9, size = 4, shape = 21) +
  theme_classic() +
  scale_shape_manual("", values = shapes[1:2]) +
  scale_fill_manual("", values = cols[1:2]) + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.position = "none")

blg.p

#ggsave("HEU_manuscript_analysis/figures/microbiome/blg_heu_nmds.pdf", device = "pdf", dpi = 300, width = 5, height = 4)

# SAF ord ####
saf.dat <- get.nmds.data(saf)

saf.p <- ggplot(saf.dat, aes(x = NMDS1, y =NMDS2, fill = site.exposure, shape = site.exposure)) +
  geom_point(alpha = 0.9, size = 4, shape = 22) +
  theme_classic() +
  scale_shape_manual("", values = shapes[5:6]) +
  scale_fill_manual("", values = cols[5:6]) + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.position = "none")

saf.p

#ggsave("HEU_manuscript_analysis/figures/microbiome/saf_heu_nmds.pdf", device = "pdf", dpi = 300, width = 5, height = 4)

# Adonis Test for Clustering ####
# determine the variance explained by HIV exposure for each site separately

physeqs <- list(cad, blg, saf)
names(physeqs) <- c("cad", "blg", "saf")


adonis.tests <- llply(as.list(physeqs), function(i){
  df <- data.frame(sample_data(i))
  dist <- phyloseq::distance(i, "bray")
  adonis(dist ~ Exposure, data = df)
})

names(adonis.tests) <- names(physeqs)

adonis.tests

# CAD R2 = 0.03527, p = 0.043
# BLG R2 = 0.03824, p = 0.579
# SAF R2 = 0.0724, p = 0.311

# END ####

