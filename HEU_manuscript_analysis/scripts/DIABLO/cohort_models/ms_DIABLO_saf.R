###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Prepare data matrices for DIABLO:South Africa
###################

# In this script, we will perform DIABLO intgration for Luminex-cytokine data and OTU data for the SOUTH AFRICAN cohort.

# load packages
library(mixOmics)
library(ggplot2)
library(plyr)
library(dplyr)

# load data
spls.dat <- readRDS("HEU_manuscript_analysis/Rdata/R_export/spls_data.rds")

# prepare data for DIABLO ####

X.otu <- spls.dat$SAF$Xotu
Y.lmx <- log10(spls.dat$SAF$Ylmx)
meta <- spls.dat$SAF$meta
Y = meta$Group

dat <- list(otu = X.otu, lmx = Y.lmx)

# check matching rownames
length(which(rownames(X.otu) ==  rownames(Y.lmx))) # all good

# set design ####

# full design with zeros diag
full.design.0 = matrix(1, ncol = length(dat), nrow = length(dat), 
                       dimnames = list(names(dat), names(dat)))
diag(full.design.0) = 0

# tune diablo ####
# set n to test
# set n to test
test.keepX = list (otu = c(5:9, seq(10, 16, 2)),
                   lmx = c(5:9, seq(10, 16, 2)))
# run tune.0
tune.saf.0 <- tune.block.splsda(X = dat, 
                                Y = Y,
                                ncomp = 2,
                                test.keepX = test.keepX,
                                design = full.design.0,
                                validation = "loo",
                                dist = "max.dist") # use loo instead of Mfold due to low n per group

#choice.keepx.0 <- tune.saf.0$choice.keepX # this varies slightly.  Here we select the number passed to the final analysis.
choice.keepx.0 <- list("otu" = c(5, 10), "lmx" = c(6, 12))

res.final <- block.splsda(X = dat,
                          Y = Y,
                          design = full.design.0,
                          keepX = choice.keepx.0,
                          scale = TRUE,
                          near.zero.var = TRUE)

# Plot Results ####
# Ord

plotIndiv(res.final, legend = TRUE, ind.names = FALSE, ellipse = TRUE) # only Comp1 is discriminatory

# Loadings ####
# comp 1
plotLoadings(res.final, comp = 1, contrib = "max", title = "comp.1")

# comp 2
plotLoadings(res.final, comp = 2, contrib = "max", title = "comp.2")

# Model error ####

error <- perf(res.final, validation = "loo")
error$error.rate.per.class 
error$error.rate # Error rate of 0.8 over two components for OTU

# save data ####

# We will save the model to compare the error from integration to sPLS-DA.
#saveRDS(res.final, "HEU_manuscript_analysis/Rdata/R_export/saf_sgccda_res.rds")

