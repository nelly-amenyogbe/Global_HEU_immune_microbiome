######################
# Nelly Amenyogbe
# HEU Manuscript Analysis: linear regression for flow cytometry cellular populations
###################

# load packages
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

# load data
flow <- read.csv("HEU_manuscript_analysis/Rdata/flow/gc_heu_fct_count.csv")

# Prepare data for analysis ####
# create molten/long data frame
colnames(flow)

f.m <- melt(flow, id.vars = c("Sample", "Stimulation", "Country", "Group"), measure.vars = c("Size", "gdTcells", "Tcells", "Bcells", "cDCs", "pDCs", "Granulocytes", "Monocytes"), variable.name = "population", value.name = "count")

all.cells <- unique(f.m$population)
all.cells

# Calculate frequency of total, as frequency of cells in the size gate
f.m <- ddply(f.m, .(Sample), transform, frequency.of.total = count/count[population == "Size"]*100)


#### Linear Model ####

# regression is performed on the frequency of total.  P-values are adjusted within each cohort specifically

populations <- all.cells[2:8] # exclude size 

cell.lm <- function(site){
  
  df <- filter(f.m, Country == site)
  
  df2 <- ldply(as.list(populations), function(i){
    
    lm.dat <- filter(df, population == i)
    
    fit <- lm(frequency.of.total ~ Group, data = lm.dat)
    
    x <- summary(fit)
    R2 <- x$adj.r.squared
    
    y <- as.data.frame(x$coefficients)
    y[,1] <- round(y[,1], digits = 3)
    y[,2] <- round(y[,2], digits = 3)
    y$confint97.5 <- y[,1] + y[,2]
    y$confint2.5 <- y[,1] - y[,2]
    y$R2 <- R2
    
    y$variable <- rownames(y)
    y$population <- i
    colnames(y)[4] <- "p.value"
    y
    
  })
  
  df2$Site <- site
  df2 <- filter(df2, variable != "(Intercept)")
  df2$p.adjust <- p.adjust(df2$p.value, method = "bonferroni")
  
  return(df2)
}



cad.lm <- cell.lm("Canada")
blg.lm <- cell.lm("Belgium")
saf.lm <- cell.lm("South.Africa")

lm.res <- rbind(cad.lm, blg.lm, saf.lm)

lm.sig <- filter(lm.res, p.adjust < 0.05) # CAD gdTcells, B cells, T cells 

# Plot Results ####

# set cols
cad.cols <- c("#006d2c", "#a1d99b")
blg.cols <- c("#08519c", "#9ecae1")
saf.cols <- c("#cb181d", "#fcbba1")

cols <- c(blg.cols, cad.cols, saf.cols)

# Prepare factor levels for plotting
f.m$Group <- gsub("Control", "HUU", f.m$Group)
f.m$Country <- gsub("South.Africa", "South Africa", f.m$Country)

f.m$site.group <- paste(f.m$Country, f.m$Group)
f.m$site.group <- factor(f.m$site.group, levels = c("Belgium HUU", "Belgium HEU", "Canada HUU", "Canada HEU", "South Africa HUU", "South Africa HEU"))

# subset dataset to only data for plotting
df.plot <- filter(f.m, population %in% populations)

ggplot(df.plot, aes(x = Country, y = frequency.of.total, fill = site.group)) + 
  geom_boxplot(outlier.size = NA, color = "black") +
  geom_point(shape = 21, size = 0.8, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  scale_fill_manual("Group", values = cols) +
  facet_wrap(~population, nrow = 2, scales = "free") + 
  theme_classic() +
  labs(x = "", y = "Frequency of Total") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 12, face = "bold", vjust = 0.8),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12))

#ggsave("HEU_manuscript_analysis/figures/flow_cytometry/gcHEU_flow_populations.pdf", device = "pdf", dpi = 300, width = 7, height = 4)
