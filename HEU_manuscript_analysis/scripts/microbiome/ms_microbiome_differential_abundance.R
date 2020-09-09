###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Microbiome Differential Abundance with DESeq2
###################

# In this script, we will identify differentially-abundant OTUs between HEU and HUU children within each site separately using DESeq2

# We will also perform sPLS-DA for all significant OTUs to determine the classification accuracy of these features in a multivariate space.

# plots produced here include results for all significant DESeq2 results, and boxplots for select discriminatory taxa 

# Load packages
library(phyloseq)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(DESeq2)
library(mixOmics)

# load helper functions
source("HEU_manuscript_analysis/scripts/functions/functions_otu_transformations.R")
source("HEU_manuscript_analysis/scripts/functions/matrix_missing_values_cleanup.R")
source("HEU_manuscript_analysis/scripts/functions/function_splsda_get_features.R")

# load data
physeq <- readRDS("HEU_manuscript_analysis/Rdata/microbiome/gc_heu_phyloseq.rds")

# aesthetic settings for plotting
sample.cols <- read.csv("HEU_manuscript_analysis/Rdata/luminex/heu_pca_colours.csv")
taxa.cols <- read.csv("HEU_manuscript_analysis/Rdata/microbiome/gc_phylum_cols.csv")

# prepare data for DESeq ####
physeq <- prune_taxa(taxa_sums(physeq) > 3, physeq) # remove any OTUs with fewer than 3 counts across the dataset

# Subset to each site
cad <- subset_samples(physeq, Site == "Canada")
cad <- prune_taxa(taxa_sums(cad) > 0, cad)

blg <- subset_samples(physeq, Site == "Belgium")
blg <- prune_taxa(taxa_sums(blg) > 0, blg)

saf <- subset_samples(physeq, Site == "South Africa")
saf <- prune_taxa(taxa_sums(saf) > 0, saf)

physeqs <- list(blg, cad, saf)
names(physeqs) <- c("blg", "cad", "saf")

# Filter OTU table to only include abundant OTUs ####
# Create matrix of OTU data

otu.tabs <- llply(physeqs, function(i){
  
  tab <- t(data.frame(otu_table(i)))
  tab
  
})

# retain only OTUs present in at least 6% of samples
otu.f <- llply(otu.tabs, function(i){
  
  otu.f <- remove.low.counts(as.data.frame(i), 1, 6)
  otu.f <- otu.f$df
  otu.f
  
})

# filter physeq objects to include relevant OTUs only
blg.f <- subset_taxa(blg, OTU.name %in% colnames(otu.f$blg))
cad.f <- subset_taxa(cad, OTU.name %in% colnames(otu.f$cad))
saf.f <- subset_taxa(saf, OTU.name %in% colnames(otu.f$saf))

ps.f <- list(blg.f, cad.f, saf.f)
names(ps.f) <- names(physeqs)

# Run DESeq2 ####
# Create DESeq objects
ds.obj <- llply(ps.f, function(i){
  
  design <- ~Exposure
  
  dds <- phyloseq_to_deseq2(i, design)
  dds
  
})

# run test
ds.tests <- llply(ds.obj, function(i){
  
  ds.wald <- DESeq(i, test = "Wald", fitType = "parametric")
  ds.wald
})

# Get deseq results ####
taxtab <- data.frame(tax_table(physeq))

ds.res <- llply(ds.tests, function(i){
  
  res <- data.frame(results(i))
  res$OTU.name <- rownames(res)
  res <- join(res, taxtab, by = "OTU.name")
  res
  
})

# significant results only
sigtabs <- llply(ds.res, function(i){
  
  sigtab <- filter(i, padj < 0.05)
  sigtab
  
})

# Plot Results ####
# Genera 
phy.cols <- as.character(taxa.cols$phylum.cols)
names(phy.cols) <- taxa.cols$Phylum

# function to produce plots for differentially-abundant OTUs
# Input: DESeq2 result output, title = character title of plot

ds.genus.logchange <- function(sigtab, title){
  
  # Arrange data
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels = names(x))
  
  # Plot
  p <- ggplot(sigtab, aes(x = Genus, y = log2FoldChange, fill = Phylum)) +
    geom_point(size = 6, alpha = 0.8, color = "black", shape = 21) +
    theme_bw() +
    scale_fill_manual("Phylum", values = phy.cols) +
    theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(face = "bold"),
          axis.line = element_line(size = 0.8),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    ggtitle(title) +
    labs(x = "", y = "Log2 Fold-Change") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#08519c")
  
  p
}

# Canada
p.cad <- ds.genus.logchange(sigtabs$cad, "Canada")
p.cad

#ggsave("HEU_manuscript_analysis/figures/microbiome/cad_heu_deseq.pdf", device = "pdf", dpi = 300, width = 7, height = 7)

# Belgium
p.blg <- ds.genus.logchange(sigtabs$blg, "Belgium")
p.blg

#ggsave("HEU_manuscript_analysis/figures/microbiome/blg_heu_deseq.pdf", device = "pdf", dpi = 300, width = 6.5, height = 3.5)

# SAF
p.saf <- ds.genus.logchange(sigtabs$saf, "South Africa")
p.saf

#ggsave("HEU_manuscript_analysis/figures/microbiome/saf_heu_deseq.pdf", device = "pdf", dpi = 300, width = 6, height = 2)

# sPLS-DA for differentially abundant OTUs ####

# Prepare sPLSDA dat ####

# Make CLR data ####
# otu.f is the list of otu tables to input
# sigtabs is the list of significant otus

otu.clr <- llply(as.list(1:length(otu.f)), function(i){
  
  otu.mat <- otu.f[[i]]
  sig.res <- sigtabs[[i]]
  
  otu.tss <- TSS.transform(otu.mat)
  otu.clr <- clr(as.matrix(otu.tss))
  
  otu.sig <- otu.clr[,which(colnames(otu.clr) %in% sig.res$OTU.name)]
  otu.sig
  
  
})

names(otu.clr) <- names(otu.f)

# Add metadata ####
meta <- data.frame(sample_data(physeq))

meta.dfs <- llply(otu.clr, function(i){
  metadat <- filter(meta, X.SampleID %in% rownames(i))
  metadat <- metadat[order(metadat$X.SampleID),]
  metadat
})

# check matching rownames.  Output should equal 1.

llply(as.list(1:length(otu.clr)), function(i){
  otu <- otu.clr[[i]]
  meta <- meta.dfs[[i]]
  
  n <- length(which(rownames(otu) == meta$X.SampleID)) / length(meta$X.SampleID)
  n
}) # all good

# BLG sPLSDA ####

X.blg <- otu.clr$blg
Y.blg <- meta.dfs$blg$Exposure

dim(X.blg)

blg.tune <- tune.splsda(X = X.blg,
                        Y = Y.blg,
                        test.keepX = c(1:17),
                        ncomp = 2,
                        validation = "Mfold",
                        folds = 2)

blg.keepx <- blg.tune$choice.keepX # 15 and 15

blg.splsda <- splsda(X = X.blg,
                     Y = Y.blg,
                     keepX = blg.keepx,
                     ncomp = 2,
                     near.zero.var = TRUE)

# CAD sPLSDA ####

X.cad <- otu.clr$cad
Y.cad <- meta.dfs$cad$Exposure

dim(X.cad)

cad.tune <- tune.splsda(X = X.cad,
                        Y = Y.cad,
                        test.keepX = c(1:47),
                        ncomp = 2,
                        validation = "Mfold",
                        folds = 4)

cad.keepx <- cad.tune$choice.keepX # 40 and 1

cad.splsda <- splsda(X = X.cad,
                     Y = Y.cad,
                     keepX = cad.keepx,
                     ncomp = 2,
                     near.zero.var = TRUE)

# SAF sPLSDA ####

X.saf <- otu.clr$saf
Y.saf <- meta.dfs$saf$Exposure

dim(X.saf)

saf.tune <- tune.splsda(X = X.saf,
                        Y = Y.saf,
                        test.keepX = c(1:3),
                        ncomp = 2,
                        validation = "Mfold",
                        folds = 2)

saf.keepx <- saf.tune$choice.keepX # 1 and 2
# force keep 3
saf.keepx <- c("comp1" = 3, "comp2" = 1)

saf.splsda <- splsda(X = X.saf,
                     Y = Y.saf,
                     keepX = saf.keepx,
                     ncomp = 2,
                     near.zero.var = TRUE)
# get error for all 
spls.obj <- list(blg.splsda, cad.splsda, saf.splsda) 
names(spls.obj) <- c("blg", "cad", "saf")

# Export error data for all cohorts ####
blg.perf <- perf(blg.splsda, validation = "Mfold", folds = 2)
cad.perf <- perf(cad.splsda, validation = "Mfold", folds = 4)
saf.perf <- perf(saf.splsda, validation = "Mfold", folds = 2)

perf.res <- list(blg.perf, cad.perf, saf.perf)
names(perf.res) <- c("Belgium", "Canada", "South Africa")

perf.res.df <- ldply(as.list(1:length(perf.res)), function(i){
  
  res <- perf.res[[i]]
  cohort <- names(perf.res)[[i]]
  
  # get overall error rate
  overall <- data.frame(res$error.rate.all$BER$max.dist)
  overall$comp <- rownames(overall)
  overall$comp <- gsub(" ", ".", overall$comp)
  overall$group <- "overall"
  colnames(overall)[1] <- "value"
  
  # get balanced error rate
  class <- data.frame(res$error.rate.class$max.dist)
  class$group <- rownames(class)
  class <- melt(class, id.vars = "group")
  colnames(class)[2] <- "comp"
  class <- class[,c("value", "comp", "group")]
  
  # combine data
  df <- rbind(overall, class)
  df$cohort <- cohort
  df
  
  
})

# plot discriminatory taxa ####

# get relevant features

# features selected on comp.1

c1.feats <- llply(spls.obj, function(i){
  feats <- mc.get.features(i, "comp.1")
  feats <- filter(feats, raw.loadings != 0)
  feats
})

# plot results

# function to plot results
# loadings = data frame output for mc.get.features, no.otus = number of OTUs to plot, otu.tab = otu table to use for plotting, cohort = one if Belgium, Canada, South Africa

plot.taxa <- function(loadings, no.otus, otu.tab, cohort){
  
  # prepare loadings
  ld <- loadings
  colnames(ld)[1] <- "OTU.name"
  ld <- join(ld, taxtab, by = "OTU.name")
  ld$OTU.name.c <- paste(ld$OTU.name, ld$Genus)
  # remove any OTUs not classified at Phylum level
  ld$max.class <- substr(ld$Species, start = 1, stop = 1)
  ld <- filter(ld, max.class != "K")
  
  ld$Group <- ifelse(ld$raw.loadings > 0, "Control",
                     ifelse(ld$raw.loadings < 0, "HEU", "Unknown"))
  
  # top heu
  ld.h <- filter(ld, Group == "HEU")
  ld.h <- ld.h[c(1:no.otus),]
  
  # top control
  ld.c <- filter(ld, Group == "Control")
  ld.c <- ld.c[c(1:no.otus),]
  
  ld.all <- rbind(ld.h, ld.c)
  
  # prepare OTU data to plot
  otus <- as.data.frame(otu.tab)
  otus <- otus[,which(colnames(otus) %in% ld.all$OTU.name)]
  otus$X.SampleID <- rownames(otus)
  otu.m <- melt(otus, id.vars = "X.SampleID")
  
  otu.meta <- join(otu.m, meta, by = "X.SampleID")
  
  # add taxonomy data
  colnames(otu.meta)[2] <- "OTU.name"
  otu.meta <- join(otu.meta, ld.all, by = "OTU.name")
  otu.meta$OTU.name.c <- paste(otu.meta$OTU.name, otu.meta$Genus)
  otu.meta$OTU.name.c <- factor(otu.meta$OTU.name.c, levels = c(as.character(ld.all$OTU.name.c)))
  
  
  # set colours
  p.cols <- filter(sample.cols, site %in% cohort)
  p.cols$group <- gsub("Control", "HUU", p.cols$group)
  nm <- as.character(p.cols$group)
  group.cols <- as.character(p.cols$colour)
  names(group.cols) <- nm
  
  # change control to HUU
  otu.meta$Group <- gsub("Control", "HUU", otu.meta$Group)
  otu.meta$Group <- factor(otu.meta$Group, levels = c("HUU", "HEU"))
  
  otu.meta$Exposure <- gsub("Control", "HUU", otu.meta$Exposure)
  otu.meta$Exposure <- factor(otu.meta$Exposure, levels = c("HUU", "HEU"))
  
  # plot
  p <- ggplot(otu.meta, aes(x = OTU.name.c, y = value, fill = Exposure)) +
    geom_boxplot(outlier.size = 1, outlier.shape = 21) +
    #geom_point(size = 2, shape=21, color = "black", position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3)) +
    scale_fill_manual("", values = group.cols) +
    facet_grid(Group~., scales = "free", space = "free") +
    theme_classic() + 
    theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14)) +
    #ggtitle(paste(cohort, group)) +
    labs(x = "", y = "Normalized Counts") + 
    coord_flip()
  
  p
  
}

# Belgium
plot.taxa(c1.feats$blg, 4, otu.clr$blg, "Belgium")
#ggsave("HEU_manuscript_analysis/figures/microbiome/blg_plot_taxa_controls.pdf", device = "pdf", dpi = 300, width = 7, height = 4.5)

# CAD
# control
plot.taxa(c1.feats$cad, 4, otu.clr$cad, "Canada")
#ggsave("HEU_manuscript_analysis/figures/microbiome/cad_plot_taxa_controls.pdf", device = "pdf", dpi = 300, width = 7, height = 4.5)

# SAF
plot.taxa(c1.feats$saf, 3, otu.clr$saf, "South Africa")
#ggsave("HEU_manuscript_analysis/figures/microbiome/saf_plot_taxa_controls.pdf", device = "pdf", dpi = 300, width = 6, height = 2.8)

# save results ####

write.csv(perf.res.df, "HEU_manuscript_analysis/Rdata/R_export/heu_microbiome_splsda_error.csv")



