###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Microbiome data preparation for DIABLO analysis
###################

# In this script, we will prepare microbiome data for use with DIABLO and sPLS-DA analyses.  This includes filtering data to remove rare OTUs and normalizing the OTU data using total sum of squares followed by the centered log ratio.

# load packages
library(phyloseq)
library(plyr)
library(dplyr)

# helper functions
source("HEU_manuscript_analysis/scripts/functions/matrix_missing_values_cleanup.R")
source("HEU_manuscript_analysis/scripts/functions/functions_otu_transformations.R")

# load data
ps <- readRDS("HEU_manuscript_analysis/Rdata/microbiome/gc_heu_phyloseq.rds")

# luminex data
lmx.dfs <- readRDS("HEU_manuscript_analysis/Rdata/R_export/heu_lmx_filtered_sPLS.RDS")

lmx.all <- lmx.dfs$all

# Prepare OTU data for filtering ####

# subset microbiome data to only include subjects with luminex data
ps.s <- subset_samples(ps, InnateID %in% unique(as.character(lmx.all$Sample))) # 72 samples with overlap

# create separate phyloseq objects for every cohort
cohorts <- c("CAD", "SAF", "BLG")

gc <- ps.s

cad <- subset_samples(gc, Site == "Canada")
saf <- subset_samples(gc, Site == "South Africa")
blg <- subset_samples(gc, Site == "Belgium")

ps <- list(cad, saf, blg)
names(ps) <- cohorts

# prune each cohort to only include OTUs present in that cohort 

# OTU tables with InnateIDs as rownames
otu.tabs <- llply(ps, function(i){
  
  tab <- data.frame(t(otu_table(i)))
  rownames(tab) <- sample_data(i)$InnateID
  tab
  
})

# TSS transform
otu.tss.tabs <- llply(otu.tabs, function(i){
  
  otu.tss <- TSS.transform(i)
  otu.tss <- as.data.frame(otu.tss)
  otu.tss
  
})

# Keep OTUs present in 50% of samples

otu.rlc <- llply(otu.tss.tabs, function(i){
  
  tab <- remove.low.counts(i, 1e-5, 30)
  return(tab$df)
  
})

# save list of OTUs to keep

otus <- llply(otu.rlc, function(i){
  
  otus <- as.character(colnames(i))
  otus
  
})

llply(otu.rlc, function(i){
  
  dim(i)
  
}) # check how many OTUs were retained

overlap.otus <- Reduce(intersect, otus)
length(overlap.otus) #378 OTUs were selected for every cohort

# Prepare OTU table from all samples together
gc.otu <- as.data.frame(otu_table(gc))
gc.meta <- data.frame(sample_data(gc))
colnames(gc.otu) <- gc.meta$InnateID # change OTU sample IDs to matching identifiers for innate data

gc.otu.s <- gc.otu[which(rownames(gc.otu) %in% overlap.otus),]
gc.otu.s <- as.data.frame(t(gc.otu.s))

# add whole cohort to rest of list
otu.rlc <- c(otu.rlc, list(gc.otu.s))
names(otu.rlc)[4] <- "ALL"

# normalize OTU data ####

# CLR
otu.clr <- llply(otu.rlc, function(i){
  mat <- as.matrix(i)
  mat.clr <- clr(mat)
  mat.clr
})

# save relevant files ####
#saveRDS(otu.clr, "HEU_manuscript_analysis/Rdata/R_export/heu_otu_clr.rds")
