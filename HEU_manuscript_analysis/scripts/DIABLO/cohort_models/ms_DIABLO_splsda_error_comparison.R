#####################
# Nelly Amenyogbe
# HEU Manuscript analysis: comparison of error rates for sPLS-DA and DIABLO models for luminex cytokines and OTU microbiome data
######################

# In this script, we will load the error rate results generated for sPLS-DA and DIABLO integration results for OTU and luminex cytokines data, for each cohort separately.  We will then construct a figure comparing these directly.

# load packages
library(mixOmics)
library(plyr)
library(dplyr)
library(reshape2)

# helper function
source("HEU_manuscript_analysis/scripts/functions/function_diablo_error.R")

# load data
group.cols <- read.csv("HEU_manuscript_analysis/Rdata/luminex/heu_pca_colours.csv") # for graphing aesthetics

# DIABLO results
blg.int <- readRDS("HEU_manuscript_analysis/Rdata/R_export/blg_sgccda_res.rds")
cad.int <- readRDS("HEU_manuscript_analysis/Rdata/R_export/cad_sgccda_res.rds")
saf.int <- readRDS("HEU_manuscript_analysis/Rdata/R_export/saf_sgccda_res.rds")

# OTU: splsda error
otu.splsda <- read.csv("HEU_manuscript_analysis/Rdata/R_export/heu_microbiome_splsda_error.csv")
otu.splsda$model <- "sPLS-DA"
otu.splsda$block <- "OTU"

lmx.splsda <- read.csv("HEU_manuscript_analysis/Rdata/R_export/heu_lmx_perf.csv")
lmx.splsda$model <- "sPLS-DA"
lmx.splsda$block <- "LMX"

# get DIABLO error results ####
# BLG
error.blg <- perf(blg.int, validation = "Mfold", folds = 2)
# CAD
error.cad <- perf(cad.int, validation = "Mfold", folds = 3)
# SAF
error.saf <- perf(saf.int, validation = "loo")

# Get error dat
g.cols <- group.cols

#BLG
blg.p <- diablo.plot.error(error.blg, g.cols, "Belgium")
blg.dat <- blg.p$plot.data

#CAD
cad.p <- diablo.plot.error(error.cad, g.cols, "Canada")
cad.dat <- cad.p$plot.data

#SAF
saf.p <- diablo.plot.error(error.saf, g.cols, "South Africa")
saf.dat <- saf.p$plot.data

# plot combined error ####
# prepare data
diablo.dat <- rbind(blg.dat, cad.dat, saf.dat)

otu.add <- otu.splsda[,c("group", "comp", "value", "block", "cohort", "model")]
lmx.add <- lmx.splsda[,c("group", "comp", "value", "block", "cohort", "model")]

all.dat <- rbind(diablo.dat, otu.add, lmx.add)

# set aesthetics
all.dat$block <- gsub("LMX", "Cytokines", all.dat$block)
all.dat$block <- gsub("OTU", "Microbiome", all.dat$block)

# plot function ###
# saves time for cohort-specific plots
# Input: ste = one of Belgium, Canada, 'South Africa'
plot.err <- function(ste){
  
  dat <- filter(all.dat, cohort == ste)
  
  dat$model <- gsub("Diablo", "DIABLO", dat$model)
  
  p.cols <- filter(group.cols, site == ste)
  cols <- c(as.character(p.cols$colour), "#969696")
  names(cols) <- c(as.character(p.cols$group), "overall")
  
  dat$model <- factor(dat$model, levels = c("sPLS-DA", "DIABLO"))
  dat <- filter(dat, comp == "comp.1")
  
  p <- ggplot(dat, aes(x = model, y = value, color=group, fill= group, group=group)) + 
    geom_point(shape=21, size = 3, color="black") +
    geom_line(size=0.8, aes(linetype=group)) +
    facet_grid(~block) +
    theme_bw() +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_linetype_manual(values = c(1, 1, 2)) +
    labs(x="", y="Maximum Error Rate") +
    ylim(c(0,1)) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = NA),
          strip.text = element_text(vjust = 0, size=14), 
          legend.key.width = unit(3, "line"),
          #legend.key.height = unit(0.2, "cm"),
          legend.title = element_blank(),
          axis.text = element_text(size=14),
          legend.text = element_text(size=14),
          axis.title.y = element_text(size=14))
  
  p
  
  
}

# blg
blg.plot <- plot.err("Belgium")
blg.plot
#ggsave("HEU_manuscript_analysis/figures/DIABLO/blg_err_clr.pdf", device="pdf", dpi=300, width=6, height=3)

# cad
cad.plot <- plot.err("Canada")
cad.plot
#ggsave("HEU_manuscript_analysis/figures/DIABLO/cad_err_clr.pdf", device="pdf", dpi=300, width=6, height=3)

# saf
saf.plot <- plot.err("South Africa")
saf.plot
#ggsave("HEU_manuscript_analysis/figures/DIABLO/saf_err_clr.pdf", device="pdf", dpi=300, width=6, height=3)

# END ####