###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Prepare data matrices for DIABLO
###################

# In this script, we will prepare microbiome and luminex data for use with DIABLO and sPLS-DA analyses.  This includes ensuring that the row names for both data sets match, and that the metadata also matches the output files.

# load packages
library(plyr)
library(dplyr)
library(reshape2)
library(missForest)

# helper functions
source("HEU_manuscript_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# load data ####
# microbiome
clr.dat <- readRDS("HEU_manuscript_analysis/Rdata/R_export/heu_otu_clr.rds")

# luminex
lmx <- readRDS("HEU_manuscript_analysis/Rdata/R_export/heu_lmx_filtered_sPLS.RDS")

# metadata 
meta <- read.csv("HEU_manuscript_analysis/Rdata/metadata/gc_heu_metadata.csv")

# taxonomy
taxonomy <- read.csv("HEU_manuscript_analysis/Rdata/microbiome/gc_heu_otu_taxonomy.csv")
colnames(taxonomy)[1] <-"OTU.num"

# Prepare luminex matrices ####
# make matrices
lmx.mats <- llply(lmx, function(i){
  mat <- dcast(i, Sample ~ cyto.stim, value.var = "final.concentration")
  rownames(mat) <- mat$Sample
  mat <- mat[,-1]
})

# remove cytokines with >15% missing values
lmx.beads <- llply(lmx.mats, function(i){
  
  mat <- beads.rm(i, 15)
  
})

lmx.beads # BLG IL-23, all IL-23.

# remove these
lmx.bd.rm <- llply(lmx.mats, function(i){
  
  mat <- remove.beads(i, 15)
  mat
  
})

# remove subjects with >15% missing data
lmx.sub <- llply(lmx.bd.rm, function(i){
  mat <- which.sub.rm(i, 15)
})

lmx.sub #BLG1033, CAD1016, SAF1023/24/27/31/33

# remove subjects
lmx.mats <- llply(lmx.bd.rm, function(i){
  mat <- subj.rm(i, 15)
  mat
})

# prepare OTU data ####
# prepare OTU tables with correct taxonomy
otu.clr.tax <- llply(clr.dat, function(i){
  tax.f <- filter(taxonomy, OTU.num %in% colnames(i))
  tax.f <- tax.f[order(tax.f$OTU.num),]
  tax.f$OTU.names <- paste(tax.f$OTU.num, tax.f$Genus)
  colnames(i) <- tax.f$OTU.names
  return(i)
})

# Prepare sPLS matrices ####

names(otu.clr.tax)
names(lmx.mats)

# ensure cohorts appear in same order
lmx.mats.reordered <- lmx.mats[c(2, 3, 1, 4)]
names(lmx.mats.reordered)

# prepare metadata ####

# create metadata tables for each cohort

meta.dfs <- llply(lmx.mats.reordered, function(i){
  
  metadata <- meta
  
  meta.s <- filter(metadata, SUBJECT_ID %in% rownames(i))
  rownames(meta.s) <- meta.s$SUBJECT_ID
  meta.s <- meta.s[order(rownames(meta.s)),]
  meta.s
  
})

# compile data
spls.data <- llply(as.list(1:length(lmx.mats)), function(i){
  
  X <- otu.clr.tax[[i]]
  Y <- lmx.mats.reordered[[i]]
  meta <- meta.dfs[[i]]
  
  # get overlapping subjects
  subs <- Reduce(intersect, list(rownames(X), rownames(Y)))
  
  # Arrange matrices to have the same subjects and in the same order
  X <- X[rownames(X) %in% subs,]
  X <- X[order(rownames(X)),]
  
  Y <- Y[rownames(Y) %in% subs,]
  Y <- Y[order(rownames(Y)),]
  Y <- missForest(Y)
  Y <- Y$ximp
  
  meta <- meta[rownames(meta) %in% subs,]
  meta <- meta[order(rownames(meta)),]
  
  to.return <- list(X, Y, meta)
  names(to.return) <- c("Xotu", "Ylmx", "meta")
  return(to.return)
  
})
names(spls.data) <- names(otu.clr.tax)

check <- clr.dat$CAD

# check matching rownames: X/Y
llply(as.list(1:length(spls.data)), function(i){
  dat <- spls.data[[i]]
  X <- dat$Xotu
  Y <- dat$Ylmx
  
  test <- length(which(rownames(X) == rownames(Y))) / nrow(X)
  test
})

# check matching rownames: X and meta
llply(as.list(1:length(spls.data)), function(i){
  dat <- spls.data[[i]]
  X <- dat$Xotu
  meta <- dat$meta
  
  test <- length(which(rownames(X) == rownames(meta))) / nrow(X)
  test
})

# save data ####
#saveRDS(spls.data, "HEU_manuscript_analysis/Rdata/R_export/spls_data.rds")



