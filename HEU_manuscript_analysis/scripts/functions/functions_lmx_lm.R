#################
# Nelly Amenyogbe
# HEU Manuscript Analysis: Functions for linear regression
#################

# These functions are written as a companion to luminex-cytokine regression analysis (see script "ms_luminex_lm.R") 

# heu.lm.adjusted ####
# HIV exposure on final concentration, adjusted by unstimulated value

# Input: cohort = one of Belgium, Canada, South Africa. stim = stimulus of interest, dat = lmx dataset

# Output:  results table for every stimulus-cytokine combination by cohort.  Results only saved for HIV exposure.

heu.lm.adjusted <- function(cohort, stim, dat){
  
  df <- filter(dat, Site == cohort, Stim == stim)
  
  df2 <- ldply(as.list(cytokines), function(i){
    
    d <- filter(df, Bead.Name == i)
    
    # run lm
    fit <- lm(log10(final.concentration) ~ Group + log10(unstim.value), data = d)
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
    
    # Filter results to only inculde HIV exposure
    y <- filter(y, variable == "GroupHEU")
    y
    
  })
  
  res <- df2
  res$Site <- cohort
  res$Stim <- stim
  
  res$Fold.Change <- res$Estimate
  res$Fold.Change <- ifelse(res$Fold.Change > 0, "HEU", "Control")
  
  return(res)
  
}

# heu.lm.adjusted.rec ####

# Function to perform regression of effect of HIV exposure on cytokine, adjusted by unstimulated value and recruitment day

# Input: cohort = one of Belgium, Canada, South Africa. stim = stimulus of interest (as indicated in data frame;), dat = lmx dataset to input

# Output:  results table for every stimulus-cytokine combination by cohort.  Results only saved for HIV exposure.

heu.lm.adjusted.rec <- function(cohort, stim, dat){
  
  df <- filter(dat, Site == cohort, Stim == stim)
  
  df2 <- ldply(as.list(cytokines), function(i){
    
    d <- filter(df, Bead.Name == i)
    
    # run lm
    fit <- lm(log10(final.concentration) ~ Group + recruit.day + log10(unstim.value), data = d)
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
    
    # Filter results to only include HIV exposure
    y <- filter(y, variable == "GroupHEU")
    y
    
  })
  
  res <- df2
  res$Site <- cohort
  res$Stim <- stim
  
  res$Fold.Change <- res$Estimate
  res$Fold.Change <- ifelse(res$Fold.Change > 0, "HEU", "Control")
  
  return(res)
  
}

# Unstimlulated function ####

# heu.lm.unstim.rec ####

# Function to perform regression of effect of HIV exposure on cytokine, adjusted by recruitment day

# Input: cohort = one of Belgium, Canada, South Africa. cyto.ex = cytokines to EXCLUDE from this analysis, dat = lmx dataset to input

# Output:  results table for every stimulus-cytokine combination by cohort.  Results only saved for HIV exposure.

heu.lm.unstim.rec <- function(cohort, cyto.ex, dat){
  
  
  df <- filter(dat, Site == cohort, Stim == "Unstim")
  df$cyto.stim <- paste0(df$Bead.Name, df$Stim)
  
  # exclude necessary cytokines
  df.ex <- df[-which(df$cyto.stim %in% cyto.ex),]
  cytos <- unique(as.character(df.ex$Bead.Name))
  
  df2 <- ldply(as.list(cytos), function(i){
    
    d <- filter(df, Bead.Name %in% i)
    
    # run lm
    fit <- lm(log10(final.concentration) ~ Group + recruit.day, data = d)
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
    y
    
  })
  
  res <- filter(df2, variable == "GroupHEU")
  res$Site <- cohort
  res$Stim <- "Unstim"
  
  res$Fold.Change <- res$Estimate
  res$Fold.Change <- ifelse(res$Fold.Change > 0, "HEU", "Control")
  
  return(res)
  
}

# heu.lm.unstim ####

# Function to perform regression of effect of HIV exposure on cytokine

# Input: cohort = one of Belgium, Canada, South Africa. cyto.ex = cytokines to EXCLUDE from this analysis, dat = lmx dataset to input

# Output:  results table for every stimulus-cytokine combination by cohort.  Results only saved for HIV exposure.

heu.lm.unstim <- function(cohort, cyto.ex, dat){
  
  
  df <- filter(dat, Site == cohort, Stim == "Unstim")
  df$cyto.stim <- paste0(df$Bead.Name, df$Stim)
  
  # exclude necessary cytokines
  df.ex <- df[-which(df$cyto.stim %in% cyto.ex),]
  cytos <- unique(as.character(df.ex$Bead.Name))
  
  df2 <- ldply(as.list(cytos), function(i){
    
    d <- filter(df, Bead.Name %in% i)
    
    # run lm
    fit <- lm(log10(final.concentration) ~ Group, data = d)
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
    y
    
  })
  
  res <- filter(df2, variable == "GroupHEU")
  res$Site <- cohort
  res$Stim <- "Unstim"
  
  res$Fold.Change <- res$Estimate
  res$Fold.Change <- ifelse(res$Fold.Change > 0, "HEU", "Control")
  
  return(res)
  
}

# Plotting lm results ####

# plot.heu ####

# This function will plot luminex cytokine results from the regression, organized by cytokines higher in either HEU or HUU children. This function uses the table summarizing significant results only.

# Input:  cohort = one of Belgium, Canada, South Africa

# Output:  ggplot boxplot overlaid with data points

plot.heu <- function(cohort, group){
  
  
  # get significant data
  sig <- filter(adj.sigtab, Site == cohort, Fold.Change == group)
  df <- filter(lmx.meta, Site == cohort, stim.cyto %in% sig$stim.cyto)
  
  # change control to HUU
  df$Group <- gsub("Control", "HUU", df$Group)
  df$Group <- factor(df$Group, levels = c("HUU", "HEU"))
  
  p.cols <- cols
  p.cols$group <- gsub("Control", "HUU", p.cols$group)
  p.cols$site.heu <- paste(p.cols$site, p.cols$group)
  
  # set aesthetics
  group.cols <- filter(p.cols, site == cohort)
  g.cols <- as.character(group.cols$colour)
  names(g.cols) <- group.cols$group
  
  
  # plot
  p <- ggplot(df, aes(x = Bead.Name, y = log10(1+final.concentration), fill = Group)) +
    geom_boxplot(color = "black", outlier.size = NA, aes(fill = Group)) +
    geom_point(alpha = 0.8, shape = 21, color = "black", size = 1, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), aes(fill = Group)) +
    facet_grid(~Stim, space = "free", scales = "free") +
    scale_color_manual(values = g.cols) +
    scale_fill_manual(values = g.cols) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0, size = 14),
          axis.text.y = element_text(size = 12),
          axis.line = element_line(size = 0.8),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12)) +
    labs(x = "", y = "Log10 Concentration") #+
  #ggtitle(cohort)
  
  return(p)
  
}
