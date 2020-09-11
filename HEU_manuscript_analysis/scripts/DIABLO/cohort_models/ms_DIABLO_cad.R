###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Prepare data matrices for DIABLO:CANADA
###################

# In this script, we will perform DIABLO intgration for Luminex-cytokine data and OTU data for the CANADIAN cohort.

# load packages
library(mixOmics)
library(ggplot2)
library(plyr)
library(dplyr)

spls.dat <- readRDS("HEU_manuscript_analysis/Rdata/R_export/spls_data.rds")

# prepare data for DIABLO ####

X.otu <- spls.dat$CAD$Xotu
Y.lmx <- log10(spls.dat$CAD$Ylmx)
meta <- spls.dat$CAD$meta
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
test.keepX = list (otu = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   lmx = c(5:9, seq(10, 18, 2), seq(20,30,5)))

# run tune.0
tune.cad.0 <- tune.block.splsda(X = dat, 
                                Y = Y,
                                ncomp = 2,
                                test.keepX = test.keepX,
                                design = full.design.0,
                                validation = "Mfold",
                                folds = 3,
                                nrepeat = 5,
                                dist = "max.dist")

#choice.keepx.0 <- tune.cad.0$choice.keepX # this varies slightly.  Here we select the number passed to the final analysis.
choice.keepx.0 <- list("otu" = c(5, 5), "lmx" = c(25, 18))

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

error <- perf(res.final, validation = "Mfold", folds = 3)
error$error.rate.per.class 
error$error.rate # Error rate of 0.44 over two components for OTU
 
# save data ####

# We will save the model to compare the error from integration to sPLS-DA.

#saveRDS(res.final, "HEU_manuscript_analysis/Rdata/R_export/cad_sgccda_res.rds")

