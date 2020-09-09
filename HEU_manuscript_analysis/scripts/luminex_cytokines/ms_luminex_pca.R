######################
# Nelly Amenoyogbe
# 09-July-2020
# Luminex: Principal Component Analysis for all ligands together
######################

# In this sript, we prepare PCAs of all luminex data:
# 1.  Coloured by Stim
# 2.  Coloured by HEU/HUU, highlighting one stim at a time with unstim
# 3. Perform the adonis test to determine variance explained by stimulus, site, and cohort of origin

# load data
lmx <- read.csv("HEU_manuscript_analysis/Rdata/luminex/gc_heu_luminex.csv")
cols <- read.csv("HEU_manuscript_analysis/Rdata/luminex/heu_pca_colours.csv")
cols$shape.fill <- c(21, 21, 24, 24, 22, 22)

# load packages
library(plyr)
library(dplyr)
library(reshape2)
library(devtools)
library(ggbiplot)
library(missForest)
library(vegan)

# load functions
source("HEU_manuscript_analysis/scripts/functions/multiplot.R") 
source("HEU_manuscript_analysis/scripts/functions/function_pca_labs.R")
source("HEU_manuscript_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# Change control to HUU ####
lmx$Group <- gsub("Control", "HUU", lmx$Group)
cols$group <- gsub("Control", "HUU", cols$group)
cols$site.heu <- paste(cols$site, cols$group)

# Make a wide dataset ####
# require(reshape2)
lmx.w <- dcast(lmx, Sample + Stim + Group + Site ~ Bead.Name, value.var = "final.concentration")

#### remove cytokines from each DF where >25% of the values are "NA"
# Prepare matrix
lmx.mat <- lmx.w
rownames(lmx.mat) <- paste0(lmx.mat$Sample, ".", lmx.mat$Stim)
lmx.w$row.names <- paste0(lmx.w$Sample, ".", lmx.w$Stim)

# remove columns with metadata
lmx.mat <- lmx.mat[,-c(1:4)]

# filter data to remove features with missing values

beads.rm(lmx.mat, percent.cutoff = 30) # no analytes to remove

# Remove subjects with > 30% NA values ####

which.s.rm <- which.sub.rm(lmx.mat, 30) # 6 subject-stimuli

# remove these from the data matrix
lmx.mat <- subj.rm(lmx.mat, 30)

# replace NA values ####
# no NA values are allowed for PCA.  Thus, I replace all NA values by imputing them with missForest.

lmx.final <- missForest(lmx.mat)
lmx.final <- lmx.final$ximp

# Make PCA of all stims #####

# create metadata file
meta <- lmx.w[,c(1:4, 19)] # we need a separate file that tells us how to colour our PCA
meta <- meta[which(meta$row.names %in% rownames(lmx.final)),]
meta$Stim <- factor(meta$Stim, levels = c("MDP", "PGN", "PAM", "LPS", "pIC", "R848", "Unstim"))

# run PCA
lmx.pca <- prcomp(log2(1 + lmx.final), scale. = TRUE)

# extract data
pca.dat <- data.frame(lmx.pca$x)
pca.dat$row.names <- rownames(pca.dat)
pca.dat <- join(pca.dat, meta, by = "row.names")
pca.dat$Stim <- factor(pca.dat$Stim, levels = c("Unstim", "PGN", "PAM", "pIC", "LPS", "R848"))

# set aesthetics for plotting
stim.cols <- c("#525252",'#beaed4','#88419d','#ffd92f','#386cb0','#f0027f','#252525')

# get axis labels
labs <- pca.labs(lmx.pca)

p.all.stims <- ggplot(pca.dat, aes(x = PC1, y = PC2, fill = Stim)) +
  geom_point(size = 3, shape = 21, color = "white") +
  scale_fill_manual("Stim", values = stim.cols) +
  theme_classic() +
  labs(x = labs$pc1, y = labs$pc2) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18),
        axis.line.x = element_line(colour = 'black', size=0.8, linetype='solid'),
        axis.line.y = element_line(colour = "black", size = 0.8, linetype = "solid"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold")) +
  ggtitle("All Ligands")

p.all.stims

#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_stims.pdf", device = "pdf", dpi = 300, width = 4, height = 4)

# Make PCA coloured by stimulus ####
# create metadata for stim-specific attributes
meta$site.heu <- paste(meta$Site, meta$Group)
meta$site.heu.stim <- paste(meta$Site, meta$Group, meta$Stim)

meta.u <- unique(meta[,c("Stim", "Group", "Site", "site.heu", "site.heu.stim")])

meta.u <- join(meta.u, cols, by = "site.heu")
meta.u$colour <- as.character(meta.u$colour)

# PCA by Stim ####

pca.cols <- function(stim){
  
  meta.u$col <- meta.u$colour
  
  meta.u$col <- ifelse(meta.u$Stim != stim, "#f0f0f0", meta.u$col)
  meta.u$col <- ifelse(meta.u$Stim == "Unstim" & meta.u$Group == "Control", "#525252", meta.u$col)
  meta.u$col <- ifelse(meta.u$Stim == "Unstim" & meta.u$Group == "HEU", "#969696", meta.u$col)
  
  # order with unstim, and then stim of interest last in the data frame
  stim.us <- filter(meta.u, Stim %in% c(stim, "Unstim"))
  rest <- filter(meta.u, Stim != stim, Stim != "Unstim")
  
  meta.reordered <- rbind(stim.us, rest)
  
  stims <- data.frame(unique(meta.reordered$Stim),"#f0f0f0")
  colnames(stims) <- c("site.heu.stim", "col")
  stims$shape.fill <- NA
  
  # combine with stimulus ellipse colours
  col <- meta.reordered[,c("site.heu.stim","col", "shape.fill")]
  col <- rbind(col, stims)
  col$col <- ifelse(col$site.heu.stim == stim, "#cb181d", col$col)
  col$col <- ifelse(col$site.heu.stim == "Unstim", "#969696", col$col)
  
  # Plot ###
  pca.dat <- data.frame(lmx.pca$x)
  pca.dat$row.names <- rownames(pca.dat)
  pca.dat <- join(pca.dat, meta, by = "row.names")
  
  # set factor levels
  pca.dat$Stim <- factor(pca.dat$Stim, levels = c(col$site.heu.stim[37:42]))
  pca.dat$site.heu.stim <- factor(pca.dat$site.heu.stim, levels = c(col$site.heu.stim[1:36]))
  
  stim.dat <- filter(pca.dat, Stim %in% c(stim, "Unstim"))
  rest.dat <- filter(pca.dat, Stim != stim, Stim != "Unstim")
  dat.re <- rbind(rest.dat, stim.dat)
  
  # Plot
  p <- ggplot(dat.re, aes(x = PC1, y = PC2, fill = site.heu.stim, shape = site.heu.stim)) +
    #stat_ellipse(aes(group = Stim, colour = Stim)) +
    geom_point(size = 3, color = "white") +
    scale_shape_manual("Group", values = c(col$shape.fill[1:36])) +
    scale_fill_manual("Group", values = c(col$col[1:36])) +
    scale_color_manual("Stim", values = c(col$col[37:42])) + 
    theme_classic() +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 15),
          axis.line.x = element_line(colour = 'black', size=0.8, linetype='solid'),
          axis.line.y = element_line(colour = "black", size = 0.8, linetype = "solid"),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14)) +
    ggtitle(stim) +
    labs(x = "", y = "")
  
  return(p)
  
  
} # creates a PCA plot coloured by Stim

# generate list of all stims
stims <- as.character(unique(meta.u$Stim))

plot.list <- llply(as.list(stims), function(i){
  
  p <- pca.cols(i)
  
})
names(plot.list) <- stims

# plot PCAs ####
# PGN
plot.list$PGN
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_pgn.pdf", device = "pdf", dpi = 300, width = 4.8, height = 4)

# PAM
plot.list$PAM 
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_pam.pdf", device = "pdf", dpi = 300, width = 4.8, height = 4)

# pIC
plot.list$pIC
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_pic.pdf", device = "pdf", dpi = 300, width = 4.8, height = 4)

# LPS
plot.list$LPS
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_lps.pdf", device = "pdf", dpi = 300, width = 4.5, height = 4)

# R848
plot.list$R848
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_r848.pdf", device = "pdf", dpi = 300, width = 4.8, height = 4)

# make individual PCAs within each stim-site ####

pca.by.stim <- function(stim, cohort){
  
  m <- filter(meta, Stim == stim, Site == cohort)
  p.d <- lmx.final[which(rownames(lmx.final) %in% m$row.names),]
  
  # ensure samples are in same order
  m <- m[order(m$row.names),]
  p.d <- p.d[order(rownames(p.d)),]
  
  # make PCA
  stim.pca <- prcomp(log2(1+p.d), scale. = TRUE)
  stim.dat <- data.frame(stim.pca$x)
  stim.dat$row.names <- rownames(stim.dat)
  stim.dat <- join(stim.dat, m, by = "row.names")
  
  # set colors
  stim.aes <- filter(cols, site == cohort)
  stim.cols <- as.character(stim.aes$colour)
  names(stim.cols) <- stim.aes$site.heu
  stim.shapes <- stim.aes$shape.fill
  names(stim.shapes) <- stim.aes$site.heu
  
  labs <- pca.labs(stim.pca)
  
  p.stim <- ggplot(stim.dat, aes(x = PC1, y = PC2, fill = site.heu, shape = site.heu)) +
    geom_point(size = 3, color = "white") +
    scale_shape_manual("Group", values = stim.shapes) +
    scale_fill_manual("Group", values = stim.cols) +
    labs(x = labs$pc1, y = labs$pc2) +
    theme_classic() +
    ggtitle(stim) +
    theme(legend.direction = "horizontal",
          legend.position = "none",
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 18),
          axis.line.x = element_line(colour = 'black', size=0.8, linetype='solid'),
          axis.line.y = element_line(colour = "black", size = 0.8, linetype = "solid"),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14, face = "bold"))
  
  return(p.stim)
  
}

cohorts <- unique(as.character(meta$Site))
stims <- unique(as.character(meta$Stim))
stims <- c("Unstim", "PGN", "PAM", "pIC", "LPS", "R848")

# CAD pcas
cad.pcas <- llply(as.list(stims), function(i){
  p <- pca.by.stim(i, "Canada")
})

names(cad.pcas) <- stims

multiplot(plotlist = cad.pcas, cols = 6)
# export as PDF, 17 x 3 inches, portrait, HEU_manuscript_analysis/figures/luminex_cytokines/cad_lmx_pca

# SAF pcas
saf.pcas <- llply(as.list(stims), function(i){
  p <- pca.by.stim(i, "South Africa")
})

names(saf.pcas) <- stims

multiplot(plotlist = saf.pcas, cols = 6)
# export as PDF, 17 x 3 inches 17 x 3 inches, portrait, HEU_manuscript_analysis/figures/luminex_cytokines/saf_lmx_pca

# BLG pcas
blg.pcas <- llply(as.list(stims), function(i){
  p <- pca.by.stim(i, "Belgium")
})

names(blg.pcas) <- stims

multiplot(plotlist = blg.pcas, cols = 6)
# export as PDF, 17 x 3 inches, portrait, HEU_manuscript_analysis/figures/luminex_cytokines/blg_lmx_pca

# lps data only for Figure 1
lps.pcas <- llply(as.list(cohorts), function(i){
  p <- pca.by.stim("LPS", i) + ggtitle("")
})
names(lps.pcas) <- cohorts

multiplot(plotlist = lps.pcas, cols = 3)
# export 9 by 3 inches, portrait, HEU_manuscript_analysis/figures/luminex_cytokines/lps_pcas

# Make legends ####
# it is not possible to generate the desired legends within the PCA plots below. Thus, there figures will incorperate these legends, which can be used with the figures above

lmx.lps <- filter(lmx, Stim == "LPS")
lmx.lps$site.heu <- paste(lmx.lps$Site, lmx.lps$Group)

legend.cols <- as.character(cols$colour)
names(legend.cols) <- as.character(cols$site.heu)


cols$unstim.cols <- c("#525252", "#969696", "#525252", "#969696", "#525252", "#969696")
unstim.cols <- as.character(cols$unstim.cols)
names(unstim.cols) <- as.character(cols$site.heu)

legend.shapes <- cols$shape
names(legend.shapes) <- as.character(cols$site.heu)

# stim cols
lmx.lps$site.heu <- factor(lmx.lps$site.heu, levels  = c(as.character(cols$site.heu)))

ggplot(lmx.lps, aes(x = Group, y = final.concentration, fill = site.heu, shape = site.heu)) +
  geom_point(size = 4) +
  theme_classic() +
  scale_shape_manual("", values = c(cols$shape.fill)) +
  scale_fill_manual("", values = legend.cols) +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12))

#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_stim_legend.pdf", device = "pdf", dpi = 300, width = 7.56, height = 4.21)

# unstim cols
ggplot(lmx.lps, aes(x = Group, y = final.concentration, fill = site.heu, shape = site.heu)) +
  geom_point(size = 4) +
  theme_classic() +
  scale_shape_manual("", values = c(cols$shape.fill)) +
  scale_fill_manual("", values = unstim.cols) +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12))

#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/heu_pca_unstim_legend.pdf", device = "pdf", dpi = 300, width = 7.56, height = 4.21)

# PCA STATS ####

# here we will perform the PERMANOVA test to assess variance explained by stimulus, site, and HIV exposure

pca.test <- log2(1+lmx.final)
test <- adonis(pca.test ~ Stim, data = meta)
test
# Stim R2 = 0.49; p < 0.001

# within stims, variance explained by site ####

ad.site <- function(stim){
  
  m <- filter(meta, Stim == stim)
  p.d <- lmx.final[which(rownames(lmx.final) %in% m$row.names),]
  
  # ensure samples are in same order
  m <- m[order(m$row.names),]
  p.d <- p.d[order(rownames(p.d)),]
  
  # run adonis
  t <- adonis(log2(1+p.d) ~ Site + Group + Site*Group, data = m)
  res <- data.frame(t$aov.tab)
  res$variable <- rownames(res)
  res$Stim <- stim
  res
  
}

stims <- unique(as.character(meta$Stim))

ad.site.res <- ldply(as.list(stims), function(i){
  
  res <- ad.site(i)
  res
  
})

# configure significant results
ad.site.res <- filter(ad.site.res, variable %in% c("Site", "Group", "Site:Group"))
colnames(ad.site.res)[6] <- "p.value"
ad.site.res$Sig <- ifelse(ad.site.res$p.value < 0.05, "Yes", "No")
ad.site.res$Sig <- factor(ad.site.res$Sig, levels = c("Yes", "No"))
ad.site.res$variable <- gsub("Group", "HIV exposure", ad.site.res$variable)
ad.site.res$Stim <- factor(ad.site.res$Stim, levels = c("Unstim", "PGN", "PAM", "pIC", "LPS", "R848"))
ad.site.res$variable <- factor(ad.site.res$variable, levels = c("Site", "HIV exposure"))

# plot R2 values for site and HIV exposure
p.ad.site <- ggplot(ad.site.res, aes(x = Stim, y = R2, color = variable, shape = Sig)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual("", values = c("#e41a1c", "#377eb8")) +
  theme(legend.position = "right",
        #legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "", y = "R2")

p.ad.site

#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/pca_ad_res_stim.pdf", device = "pdf", dpi = 300, width = 5.5, height = 2.3)

# Within stims & sites, HIV exposure ###

ad.group <- function(site){
  
  m <- filter(meta, Site == site)
  p.d <- lmx.final[which(rownames(lmx.final) %in% m$row.names),]
  
  ad.res <- ldply(as.list(stims), function(i){
    
    m <- filter(m, Stim == i)
    p.d <- p.d[which(rownames(p.d) %in% m$row.names),]
    
    # ensure samples are in same order
    m <- m[order(m$row.names),]
    p.d <- p.d[order(rownames(p.d)),]
    
    # run adonis
    t <- adonis(log2(1+p.d) ~ Group, data = m)
    res <- data.frame(t$aov.tab)
    res$variable <- rownames(res)
    res$Stim <- i
    res
    
  })
  
  ad.res$Site <- site
  ad.res <- filter(ad.res, variable == "Group")
  ad.res
  
} 

sites <- unique(as.character(meta$Site))

ad.res.all <- ldply(as.list(sites), function(i){
  
  res <- ad.group(i)
  res
  
})

# configure adonis results
colnames(ad.res.all)[6] <- "p.value"
ad.res.all$Sig <- ifelse(ad.res.all$p.value < 0.05, "Yes", "No")
ad.res.all$Sig <- factor(ad.res.all$Sig, levels = c("Yes", "No"))

ad.res.all$Stim <- factor(ad.res.all$Stim, levels = c("Unstim", "PGN", "PAM", "pIC", "LPS", "R848"))

# graph results:  HEU exposure
site.cols <- filter(cols, group == "HUU")
s.cols <- as.character(site.cols$colour) 
names(s.cols) <- site.cols$site

p.ad <- ggplot(ad.res.all, aes(x = Stim, y = R2, colour = Site, shape = Sig)) + geom_point(size = 3) +
  ylim(0,0.3) +
  theme_bw() +
  scale_color_manual(values = s.cols) +
  theme(legend.position = "right",
        #legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = "", y = "R2")

p.ad
#ggsave("HEU_manuscript_analysis/figures/luminex_cytokines/pca_ad_res.pdf", device = "pdf", dpi = 300, width = 5.5, height = 2.3)
