##################
# Nelly Amenyogbe
# 12-Sept-2018
# function: plot error for diablo results
################

# This function was written specifically to plot the error rate for a given DIABLO result for each data block separately.  It is customized for the data provided for the HEU Manuscript Analysis. This is applied in the script DIABLO/ms_DIABLO_splsda_error_comparison.R.

# Input ####
# perf.obj  = the output of the perf() function implemented by mixOmics v6, applied to a block.splsda result.  

# group.cols <- data frame supplied by read.csv("Rdata/luminex/heu_pca_colours.csv)

# cohort <- one of "Belgium", "Canada", "South Africa"

# Output ####

# error.plot: ggplot-generated line plot assessing max error rate for comp 1 and comp 2, for cytokine (LMX) and microbiome (OTU) data

# plot.data: data frame used to produce the error.plot ggplot object


diablo.plot.error <- function(perf.obj, group.cols, cohort){
  
  overall.otu <- data.frame(perf.obj$error.rate$otu)
  overall.otu$comp <- rownames(overall.otu)
  overall.otu$comp <- gsub(" ", ".", overall.otu$comp)
  overall.otu$block <- "OTU"
  overall.otu$group <- "overall"
  overall.otu <- overall.otu[,c("group", "comp", "max.dist", "block")]
  colnames(overall.otu)[2:3] <- c("variable", "value")
  
  overall.lmx <- data.frame(perf.obj$error.rate$lmx)
  overall.lmx$comp <- rownames(overall.lmx)
  overall.lmx$comp <- gsub(" ", ".", overall.lmx$comp)
  overall.lmx$block <- "LMX"
  overall.lmx$group <- "overall"
  overall.lmx <- overall.lmx[,c("group", "comp", "max.dist", "block")]
  colnames(overall.lmx)[2:3] <- c("variable", "value")
  
  class.otu <- data.frame(perf.obj$error.rate.per.class$otu$max.dist)
  class.otu <- class.otu[,c(1:2)]
  class.otu$group <- rownames(class.otu)
  colnames(class.otu)[1:2] <- c("comp.1", "comp.2")
  otu.m <- melt(class.otu, id.vars = "group")
  otu.m$block <- "OTU"
  
  class.lmx <- data.frame(perf.obj$error.rate.per.class$lmx$max.dist)
  class.lmx <- class.lmx[,c(1:2)]
  class.lmx$group <- rownames(class.lmx)
  colnames(class.lmx)[1:2] <- c("comp.1", "comp.2")
  lmx.m <- melt(class.lmx, id.vars = "group")
  lmx.m$block <- "LMX"
  
  error.dat <- rbind(overall.otu, overall.lmx, otu.m, lmx.m)
  error.dat$group <- gsub("UE", "Control", error.dat$group)
  error.dat$cohort <- cohort
  error.dat$model <- "Diablo"
  colnames(error.dat)[2] <- "comp"
  
  p.cols <- filter(group.cols, site == cohort)
  cols <- c(as.character(p.cols$colour), "#969696")
  names(cols) <- c(as.character(p.cols$group), "overall")
  
  
  p.error <- ggplot(error.dat, aes(x = comp, y = value, group = group, colour = group)) +
    geom_line(size = 0.8, aes(linetype = group)) +
    theme_bw() +
    facet_wrap(~block) +
    scale_color_manual("Group", values = cols) + 
    scale_linetype_manual("Group", values = c(1,1,1,1,1,1,5)) +
    labs(x = "", y = "Error Rate") +
    theme(legend.position = "bottom",
          strip.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12),
          axis.line = element_line(size = 0.8),
          strip.background = element_rect(fill = "white")) +
    ylim(c(0, 1))
  
  to.return <- list(p.error, error.dat)
  names(to.return) <- c("error.plot", "plot.data")
  
  return(to.return)
  
}



