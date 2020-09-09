######################
# Nelly Amenyogbe
# HEU Manuscript Analysis: linear regression for flow cytometry intracellular cytokine expression
###################

# load packages
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

# load data
df <- read.csv("HEU_manuscript_analysis/Rdata/flow/gc_heu_fct_icc.csv")
meta <- read.csv("HEU_manuscript_analysis/Rdata/metadata/gc_heu_metadata.csv")
cols <- read.csv("HEU_manuscript_analysis/Rdata/luminex/heu_pca_colours.csv") # for graphing aesthetics

# prepare data for lm test
df.test <- df

df.test$cell.cyto <- paste0(df.test$population, "_", df.test$cytokine)
df.test$cell.cyto.sub <- paste0(df.test$cell.cyto, df.test$Sample.Name)

# add unstim.value
df.test <- ddply(df.test, .(cell.cyto.sub), transform, unstim.value = value[Stim == "Unstim"])
df.test$feat.unique <- paste0(df.test$cell.cyto, "_", df.test$Stim)

# add recrutment day
rec.add <- meta[,c("SUBJECT_ID", "recruit.day")]
colnames(rec.add)[1] <- "Sample.Name"
df.test <- join(df.test, rec.add, by = "Sample.Name")


df.stim <- filter(df.test, Stim != "Unstim")
df.us <- filter(df.test, Stim == "Unstim")

# get unique variables
vars <- unique(as.character(df.stim$feat.unique))
vars.us <- unique(as.character(df.us$feat.unique))
sites <- unique(as.character(df.test$Country))


# Run LM ####
# test 

heu.lm <- function(cohort){
  
  df <- filter(df.stim, Country == cohort)
  
  df <- filter(df, population %in% c("pDC", "cDC", "Monocytes"))
  vars <- unique(as.character(df$feat.unique))
  
  df.res <- ldply(as.list(vars), function(i){
    
    lm.dat <- filter(df, feat.unique == i)
    
    fit <- lm(value ~ Group + recruit.day + unstim.value, data = lm.dat)
    
    # save results
    x <- summary(fit)
    R2 <- x$adj.r.squared
    
    y <- as.data.frame(x$coefficients)
    y[,1] <- round(y[,1], digits = 3)
    y[,2] <- round(y[,2], digits = 3)
    y$confint97.5 <- y[,1] + y[,2]
    y$confint2.5 <- y[,1] - y[,2]
    y$R2 <- R2
    
    y$variable <- rownames(y)
    y$feat.unique <- i
    colnames(y)[4] <- "p.value"
    y
    
  })
  
  df.res$Site <- cohort
  df.res <- filter(df.res, variable != "(Intercept)", variable != "unstim.value")
  
  return(df.res)
  
}

cad.lm <- heu.lm("Canada")
blg.lm <- heu.lm("Belgium")
saf.lm <- heu.lm("South Africa")

# unstim data ####

heu.lm.us <- function(cohort){
  
  df <- filter(df.us, Country == cohort)
  df <- filter(df, population %in% c("cDC", "pDC", "Monocytes"))
  vars.us <- unique(as.character(df$feat.unique))
  
  df.res <- ldply(as.list(vars.us), function(i){
    
    lm.dat <- filter(df, feat.unique == i)
    
    fit <- lm(value ~ Group, data = lm.dat)
    
    # save results
    x <- summary(fit)
    R2 <- x$adj.r.squared
    
    y <- as.data.frame(x$coefficients)
    y[,1] <- round(y[,1], digits = 3)
    y[,2] <- round(y[,2], digits = 3)
    y$confint97.5 <- y[,1] + y[,2]
    y$confint2.5 <- y[,1] - y[,2]
    y$R2 <- R2
    
    y$variable <- rownames(y)
    y$feat.unique <- i
    colnames(y)[4] <- "p.value"
    y
    
  })
  
  df.res$Site <- cohort
  df.res <- filter(df.res, variable != "(Intercept)")
  
  return(df.res)
  
}

cad.lm.us <- heu.lm.us("Canada")
blg.lm.us <- heu.lm.us("Belgium")
saf.lm.us <- heu.lm.us("South Africa")

lm.res <- rbind(blg.lm, blg.lm.us, cad.lm, cad.lm.us, saf.lm, saf.lm.us)

# finalize lm.res ####

# add descriptors of cell populations
cell.meta <- df.test[,c("feat.unique", "Stim", "cytokine", "population")]
cell.meta <- unique(cell.meta)

lm.res <- join(lm.res, cell.meta, by = "feat.unique")
lm.res <- filter(lm.res, variable == "GroupHEU")

# adjust p ####
sites <- unique(as.character(lm.res$Site))
cells <- unique(as.character(lm.res$population))

lm.res <- ldply(as.list(sites), function(i){
  
  df <- ldply(as.list(cells), function(j){
    
    dat <- filter(lm.res, Site == i, population == j)
    dat$p.adj <- p.adjust(dat$p.value, method = "BH")
    dat
    
  })
  
})



lm.res$sig <- ifelse(lm.res$p.adj < 0.001, "***",
                     ifelse(lm.res$p.adj < 0.01, "**",
                            ifelse(lm.res$p.adj < 0.05, "*",
                                   ifelse(lm.res$p.adj < 0.10, ".", ""))))


# plot results ####
df.test$stim.cyto <- paste0(df.test$Stim,"_", df.test$cytokine)

sigtab <- filter(lm.res, p.adj < 0.05) # Canada only

# Canada
# get significant findings
cad.sig <- filter(lm.res, Site == "Canada", p.adj < 0.05)

# filter data to significant features
cad.p.sig <- filter(df.test, Country == "Canada", feat.unique %in% cad.sig$feat.unique)

ggplot(cad.p.sig, aes(x = stim.cyto, y = log10(1+value), colour = Group)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(size = 0.8, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  facet_grid(~population, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0)) +
  ggtitle("Canada:  Significant results")

# Plot IL-12 for Canada ####

# subset data for plotting
il12.p <- filter(df.test, cytokine == "IL-12", Country == "Canada", population %in% c("cDC", "Monocytes"), Stim %in% c("Unstim", "PAM", "PGN", "LPS", "R848"))


il12.p$Stim <- factor(il12.p$Stim, levels = c("Unstim", "PAM", "PGN", "LPS", "R848"))

# change Control to HUU
il12.p$Group <- gsub("Control", "HUU", il12.p$Group)
il12.p$Group <- factor(il12.p$Group, levels = c("HUU", "HEU"))

p.12 <- ggplot(il12.p, aes(x = Stim, y = value, fill = Group)) +
  geom_boxplot(outlier.size = NA, color = "black") +
  facet_wrap(~population, scales = "free") +
  geom_point(shape = 21, alpha = 0.8, size = 1, position = position_jitterdodge(dodge.width = 0.7)) +
  theme_classic() +
  scale_fill_manual("", values = c("HUU" = "#006d2c", "HEU" = "#a1d99b")) +
  theme(axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ylim(c(0,100)) +
  labs(x = "", y = "% Producers")

p.12

#ggsave("HEU_manuscript_analysis/figures/flow_cytometry/cad_heu_flow_il12.pdf", device = "pdf", dpi = 300, width = 7, height = 2.5)

# plot rest of Canada sig results ###
cad.p <- filter(cad.p.sig, Country == "Canada", population %in% c("cDC", "Monocytes"), Stim %in% c("Unstim", "PAM", "PGN", "LPS", "R848"))

cad.p$Stim <- factor(cad.p$Stim, levels = c("Unstim", "PAM", "PGN", "LPS", "R848"))

p.cads <- ggplot(filter(cad.p, cytokine != "IL-12"), aes(x = Stim, y = value, fill = Group)) +
  geom_boxplot(outlier.size = NA, color = "black") +
  facet_wrap(~population + cytokine, scales = "free") +
  geom_point(shape = 21, alpha = 0.8, size = 1, position = position_jitterdodge(dodge.width = 0.7)) +
  theme_classic() +
  scale_fill_manual(values = c("Control" = "#006d2c", "HEU" = "#a1d99b")) +
  theme(axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  #ylim(c(0,100)) +
  labs(x = "", y = "% Producers")

p.cads

#ggsave("HEU_manuscript_analysis/figures/flow_cytometry/cad_heu_flow_others.pdf", device = "pdf", dpi = 300, width = 3.9, height = 2.5)
