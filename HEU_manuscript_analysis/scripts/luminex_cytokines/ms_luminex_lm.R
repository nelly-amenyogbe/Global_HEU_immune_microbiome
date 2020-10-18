###################
# Nelly Amenyogbe
# HEU Manuscript Analysis: linear regression for luminex cytokine data
###################

# In this script, we perform linear regression to determine the effect of HIV exposure on cytokine responses to PRR stimulation and create a graph of the significant results

# load packages
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

# load functions
source("HEU_manuscript_analysis/scripts/functions/functions_lmx_lm.R")

# load data
lmx <- read.csv("HEU_manuscript_analysis/Rdata/luminex/gc_heu_luminex.csv")
meta <- read.csv("HEU_manuscript_analysis/Rdata/metadata/gc_heu_metadata.csv") # Metadata
cols <- read.csv("HEU_manuscript_analysis/Rdata/luminex/heu_pca_colours.csv")

# prepare luminex data for analysis ####

# remove zero values ####

# flag anaytes with zeroes
zero <- filter(lmx, final.concentration == 0)
unique(zero$Bead.Name) # IL-23 only

# These will be replaced with 0.01 * the minimum concentration
il23 <- filter(lmx, Bead.Name == "IL-23", final.concentration > 0)
il23.min <- min(il23$final.concentration) # 1
# change zero-values to 0.01
lmx$final.concentration <- ifelse(lmx$final.concentration == 0, 0.01, lmx$final.concentration)

# Add unstim value to data ####
# this will be used to adjust for baseline values

lmx$subject.bead <- paste(lmx$Sample, lmx$Bead.Name)

lmx <- ddply(lmx, .(subject.bead), transform, unstim.value = final.concentration[Stim == "Unstim"])

# Add metadata ####
meta <- meta[,c("SUBJECT_ID", "SEX", "GEST_AGE", "BIRTH_WT", "DELIVERY", "MOM_AGE", "WAZ", "WLZ", "HAZ", "art.ga", "recruit.day")]

colnames(meta)[1] <- "Sample"

lmx.meta <- join(lmx, meta, by = "Sample")

# Write Linear Model:  All data together ####
cytokines <- unique(as.character(lmx.meta$Bead.Name))
stims <- unique(as.character(lmx.meta$Stim))
stims <- stims[2:6] # remove unstim
cohorts <- unique(as.character(lmx.meta$Site))

# Adjusted by unstimulated value only ####
# for Belgium and South Africa

res.blg.saf <- ldply(as.list(cohorts[c(1,3)]), function(i){
  
  ldply(as.list(stims), function(j){
    heu.lm.adjusted(i, j, lmx.meta)
  })
  
})

# Adjusted by unstimulated value and recruit day ####
# For Canada only

res.cad <- ldply(as.list(stims), function(i){
  
  heu.lm.adjusted.rec("Canada", i, lmx.meta)
  
})

all.adj <- rbind(res.blg.saf, res.cad)

# Unstimualed baseline analysis ####
# linear regression for unstimulated values

# Filter unstimulated values
# we will only select unstim variables where over 70% of subjects had a measured response AND the difference between TRUE and FALSE MFI is  > 10

lmx.us <- filter(lmx, Stim == "Unstim") # unstimulated data only

# flag extrapolated data as estimated values in raw concentratin
lmx.us$cutoff <- substr(lmx.us$rawConcentration, start = 1, stop = 1)
lmx.us$cutoff <- ifelse(lmx.us$cutoff == "<" | lmx.us$Bead.Name == "IL-23" & lmx.us$final.concentration == 0.04, "FALSE", "TRUE")
lmx.us$cyto.stim <- paste0(lmx.us$Bead.Name, lmx.us$Stim)

# Determine the frequency of low MFI for all unstimulated cytokines

us.res <- llply(as.list(cohorts), function(i){
  
  df <- filter(lmx.us, Site %in% c(i))
  cytokines <- unique(as.character(df$cyto.stim))
  
  df.result <- ldply(as.list(cytokines), function(i){
    
    df.test <- filter(df, cyto.stim == i)
    freq.true <- length(which(df.test$cutoff == "TRUE")) / length(df.test$cutoff)
    freq.mfi10 <- length(which(df.test$MFI < 10)) / length(df.test$MFI)
    
    cyto.stim <- i
    res <- data.frame(cyto.stim, freq.true, freq.mfi10)
  })
  
  df.result
})

names(us.res) <- cohorts

# selet unstim data to remove
us.remove <- llply(us.res, function(i){
  
  res <- i
  res <- filter(res, freq.true <= 0.7 & freq.mfi10 > 0.3)
  rm <- as.character(res$cyto.stim)
  rm
  
})

us.remove

# Run baseline lm ####

cad.us <- heu.lm.unstim.rec("Canada", us.remove$Canada, lmx.meta)
blg.us <- heu.lm.unstim("Belgium", us.remove$Belgium, lmx.meta)
saf.us <- heu.lm.unstim("South Africa", us.remove$`South Africa`, lmx.meta)

us.adj <- rbind(blg.us, cad.us, saf.us)

# lm result summarization ####
all.lm.res <- rbind(all.adj, us.adj)

# adjust p-values
# p-values adjusted for all stim-cytokine combinations within each cohort separately via benjamini-hochberg method

all.lm.res <- ldply(as.list(cohorts), function(i){
  df <- filter(all.lm.res, Site == i)
  df$p.adj <- p.adjust(df$p.value, method = "BH")
  df
})


adj.sigtab <- filter(all.lm.res, p.adj < 0.051) # one result p = 0.05 (Belgium LPS IL-6 p = 0.05004).  Retain this and present data accordingly

# save results
# write.csv(all.res, "gc_heu_lmx_lm_res.csv") # not run

# Graph significant results ####

# filter luminex data to only include significant features
adj.sigtab$stim.cyto <- paste(adj.sigtab$Stim, adj.sigtab$Bead.Name)
lmx.meta$stim.cyto <- paste(lmx.meta$Stim, lmx.meta$Bead.Name)

unique(adj.sigtab$Site) # No significant findings for South Africa

# generate plots 

# Plot Canada #### 
# Controls
plot.heu("Canada", "Control")
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/cad_lm_c.pdf", device = "pdf", dpi = 300, height = 4.0, width = 7.21)

# HEU
plot.heu("Canada", "HEU")
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/cad_lm_h.pdf", device = "pdf", dpi = 300, height = 4, width = 3.57)

# Plot Belgium #### 

# HEU (no controls)
plot.heu("Belgium", "HEU") 
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/blg_lm_h.pdf", device = "pdf", dpi = 300, height = 4, width = 4.08)

# BLG ART stratification ####
# here, we will determine if initiating ART prior to pregnancy is associated with any differences between the HEU children

# generate dataset for BLG infants only
blg.heu <- filter(lmx.meta, Site == "Belgium", Group == "HEU")

# create art initiation variable

blg.heu$initiation <- ifelse(blg.heu$art.ga == "preconception", "preconception",
                             ifelse(blg.heu$art.ga %in% c("T1", "T2", "T3"), "during_pregnancy", blg.heu$art.ga))

# check results
table(blg.heu$art.ga, blg.heu$initiation) # there is one Belgian who initiated ART after delivery.  Remove from this analysis, as not enough group membership.

blg.heu <- filter(blg.heu, art.ga != "after.delivery")

table(blg.heu$art.ga, blg.heu$initiation)

# get list of stim cytos
blg.test <- filter(adj.sigtab, Site == "Belgium")
blg.test <- blg.test$stim.cyto

# test lm with this list

blg.res <- ldply(as.list(blg.test), function(i){
  
  d <- filter(blg.heu, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ initiation + log10(unstim.value), data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable == "initiationpreconception")
  
  y
  
  
}) 

blg.res$p.value # no significant findings

# regressions adjusted for unbalanced host factors ####
# in this section, we will perform cohort-specific regressions, adjusting for host factors that differ betweeh HEU and HUU children, to determine if any significant findings are driven by host factors other than HIV exposure

# BLG:  results adjusted for sex ####

blg.dat <- filter(lmx.meta, Site == "Belgium")

blg.res.sex <- ldply(as.list(blg.test), function(i){
  
  d <- filter(blg.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + SEX + log10(unstim.value), data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
}) 

blg.res.sex$p.value[blg.res.sex$variable == "GroupHEU"] # results remain significant when adjusting for sex

# BLG:  results adjusted for DM ####

blg.res.dm <- ldply(as.list(blg.test), function(i){
  
  d <- filter(blg.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + DELIVERY + log10(unstim.value), data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
})

blg.res.dm$p.value[blg.res.dm$variable == "GroupHEU"] # results remain significant when adjusting for delivery mode

# CAD results adjusted for maternal age ####

# get significant features to test
cad.vars <- filter(adj.sigtab, Site == "Canada")
cad.vars <- cad.vars$stim.cyto

# Get Canadian data
cad.dat <- filter(lmx.meta, Site == "Canada")

# excluding unstimulated data
cad.res.ma <- ldply(as.list(cad.vars[1:8]), function(i){
  
  d <- filter(cad.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + MOM_AGE + recruit.day + log10(unstim.value), data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
}) 

cad.res.ma$p.value[cad.res.ma$variable == "GroupHEU"] # results remain significant when adjusting for maternal age

# unstimulated data only
cad.res.ma.us <- ldply(as.list(cad.vars[9:11]), function(i){
  
  d <- filter(cad.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + MOM_AGE + recruit.day, data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
}) 

cad.res.ma.us$p.value[cad.res.ma.us$variable == "GroupHEU"] # results remain significant when adjusting for maternal age

# CAD results adjusted for WLZ ####

cad.res.wlz <- ldply(as.list(cad.vars[1:8]), function(i){
  
  d <- filter(cad.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + WLZ + recruit.day + log10(unstim.value), data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
}) 

# stimulated values
cad.res.wlz$p.value[cad.res.wlz$variable == "GroupHEU"] # all significant after adjusting

# unstimulated values

cad.res.wlz.us <- ldply(as.list(cad.vars[9:11]), function(i){
  
  d <- filter(cad.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + WLZ + recruit.day, data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
}) 

cad.res.wlz.us$p.value[cad.res.wlz.us$variable == "GroupHEU"] # all significant after adjusting

# CAD results adjusted for LAZ ####

# stimulated values
cad.res.laz <- ldply(as.list(cad.vars[1:8]), function(i){
  
  d <- filter(cad.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + LAZ + log10(unstim.value), data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
})

cad.res.laz$p.value[cad.res.laz$variable == "GroupHEU"] # all significant after adjusting

# stimulated values
cad.res.laz.us <- ldply(as.list(cad.vars[9:11]), function(i){
  
  d <- filter(cad.dat, stim.cyto == i)
  
  fit <- lm(log10(final.concentration) ~ Group + LAZ + recruit.day, data = d)
  x <- summary(fit) # get results
  r2 <- x$r.squared
  
  # Create results summary
  y <- as.data.frame(x$coefficients)
  y[,1] <- round(y[,1], digits = 3)
  y[,2] <- round(y[,2], digits = 3)
  y$confint97.5 <- y[,1] + y[,2]
  y$confint2.5 <- y[,1] - y[,2]
  y$variable <- rownames(y)
  y$R2 <- r2
  y$Bead.Name <- i
  colnames(y)[4] <- "p.value"
  
  y <- filter(y, variable != "(Intercept)", variable != "log10(unstim.value)")
  
  y
  
  
}) 

cad.res.laz.us$p.value[cad.res.laz.us$variable == "GroupHEU"] # all significant after adjusting

# conclude: demographic host factors do not account for the differences between HEU and HUU children in this study