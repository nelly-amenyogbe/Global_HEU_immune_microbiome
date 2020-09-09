################
# Nelly Amenyogbe
# 15-Feb-2019
# Funtion:  pca.labs
###############

# This fuction creates axis labels for a PCA generated using the prcomp function in base R.

# Input ####
# pca object, generated with the prcomp function

# Output ####
# list; axis labels for the first and second PCA components that can be used to label the pca plot, structured as, for example, "PC1: 39.2%" 


pca.labs <- function(pca){
  
  res <- summary(pca)$importance
  pc1 <- round(100*(res[2,1]), digits = 1)
  pc2 <- round(100*(res[2,2]), digits = 1)
  pc1lab <- paste0("PC1:", " ", pc1, "%")
  pc2lab <- paste0("PC2:", " ", pc2, "%")
  
  labs <- list(pc1lab, pc2lab)
  names(labs) <- c("pc1", "pc2")
  return(labs)
  
}