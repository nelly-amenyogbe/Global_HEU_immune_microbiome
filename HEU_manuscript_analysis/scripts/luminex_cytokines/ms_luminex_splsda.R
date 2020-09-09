########################
# Nelly Amenyogbe
# HEU mausript analysis: sPLS-DA analysis for luminex cytokine data
#######################

# In this script we will perform sPLS-DA for luminex cytokine data. Error rates from the sPLS-DA models will be compared to error rates from DIABLO integrations.

# Here we will not use the same data that was generated for DIABLO as there is no need to pre-select subjects that only have microbiome data.  Thus, prior to performing sPLS-DA, we will prepare the data matrices using the dataset with "unresponsive" cytokines identified by the fligner-kileen test removed.

# load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(mixOmics)
library(reshape2)
library(missForest)

# helper functions
source("HEU_manuscript_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# load data
lmx <- readRDS("HEU_manuscript_analysis/Rdata/R_export/heu_lmx_filtered_sPLS.RDS")

# Prepare matrices for sPLS-DA ####
# data.frames with HEU status retained
lmx.cast <- llply(lmx, function(i){
  
  df.cast <- dcast(i, Sample + Group ~ cyto.stim, value.var = "final.concentration")
  rownames(df.cast) <- df.cast$Sample
  df.cast
  
})

# data matrices
mats <- llply(lmx.cast, function(i){
  
  mat <- i
  mat <- mat[,-c(1:2)]
  mat
  
})

# flag cytokines to remove
cyto.rm <- llply(mats, function(i){
  
  rm <- beads.rm(i, 15)
  rm
  
})

cyto.rm # IL-23 from Belgium

# flag subjects to remove
sub.rm <- llply(mats, function(i){
  
  rm <- which.sub.rm(i, 15)
  rm
  
}) 
sub.rm # 1 BLG, 3 SAF

# make final matrices
mats.final <- llply(mats, function(i){
  
  # remove subjects and cytokines with missing values
  mat <- i
  mat.bead.rm <- remove.beads(mat, 15)
  mat.sub.rm <- subj.rm(mat.bead.rm, 15)
  
  # replace remaining NA values
  mat.final <- missForest(mat.sub.rm)
  mat.final <- mat.final$ximp
  mat.final <- mat.final[order(rownames(mat.final)),]
  
  return(mat.final)
})

# Create metadata

meta <- llply(as.list(1:length(mats.final)), function(i){
  
  meta <- lmx.cast[[i]]
  mat <- mats.final[[i]]
  
  meta <- meta[,c(1:2)]
  meta.f <- filter(meta, Sample %in% rownames(mat))
  meta.f <- meta.f[order(meta.f$Sample),]
  
  return(meta.f)
})
names(meta) <- names(mats.final)

# check matching rownames (all should equal 1)
check <- llply(as.list(1:length(mats.final)), function(i){
  
  meta <- meta[[i]]
  mat <- mats.final[[i]]
  
  n <- length(which(rownames(mat) == meta$Sample)) / length(meta$Sample)
  n
  
})

check # all match :)

# run sPLS-DA ####

# Canada sPLS-DA ####

X.c <- log10(mats.final$Canada)
Y.c <- meta$Canada$Group

n.test <- c(1:dim(X.c)[1])

c.tune <- tune.splsda(X = X.c,
                      Y = Y.c,
                      ncomp = 2,
                      test.keepX = n.test,
                      validation = "Mfold", 
                      folds = 4)

c.keep.x <- c.tune$choice.keepX # varies each time.  Set to last one passed to analysis
c.keep.x <-  c("comp1" = 24, "comp2" = 3)

c.splsda <- splsda(X = X.c,
                   Y = Y.c,
                   ncomp = 2,
                   keepX = c(c.keep.x))

c.perf <- perf(c.splsda, validation = "Mfold", folds = 4)
plot(c.perf)

# Belgium sPLS-DA ####

X.b <- log10(mats.final$Belgium)
Y.b <- meta$Belgium$Group

# tune
n.test <- c(1:dim(X.b)[1])

b.tune <- tune.splsda(X = X.b,
                      Y = Y.b,
                      ncomp = 2,
                      test.keepX = n.test,
                      validation = "Mfold", 
                      folds = 4)

b.keep.x <- b.tune$choice.keepX # varies each time.  Set to last one passed to analysis
b.keep.x <-  c("comp1" = 9, "comp2" = 8)


b.splsda <- splsda(X = X.b,
                   Y = Y.b,
                   ncomp = 2,
                   keepX = c(b.keep.x))

b.perf <- perf(b.splsda, validation = "Mfold", folds = 4)
plot(b.perf)

# SAF sPLS-DA ####

X.s <- log10(mats.final$`South Africa`)
Y.s <- meta$`South Africa`$Group

# tune
n.test <- c(1:dim(X.s)[1])

s.tune <- tune.splsda(X = X.s,
                      Y = Y.s,
                      ncomp = 2,
                      test.keepX = n.test,
                      validation = "Mfold", 
                      folds = 4)

s.keep.x <- s.tune$choice.keepX  # varies each time.  Set to last one passed to analysis
s.keep.x <-  c("comp1" = 13, "comp2" = 2)

s.splsda <- splsda(X = X.s,
                   Y = Y.s,
                   ncomp = 2,
                   keepX = c(s.keep.x))

s.perf <- perf(s.splsda, validation = "Mfold", folds = 4)
plot(s.perf)

# Get error ####
perf.res <- list(b.perf, c.perf, s.perf)
names(perf.res) <- c("Belgium", "Canada", "South Africa")

# combined error data from all cohorts

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

# add annotation for site and HEU status
perf.res.df$col.group <- ifelse(perf.res.df$group %in% c("Control", "HEU"), paste(as.character(perf.res.df$cohort), as.character(perf.res.df$group)), as.character(perf.res.df$group))

# save error dat ####
#write.csv(perf.res.df, "HEU_manuscript_analysis/Rdata/R_export/heu_lmx_perf.csv")

