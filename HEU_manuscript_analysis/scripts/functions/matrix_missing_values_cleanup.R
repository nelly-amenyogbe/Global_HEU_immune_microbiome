#######################
# Nelly Amenyogbe
# 07 December 2016
# Functions for identifying and removing features with missing values in matrices
#######################

#### Column remove List ####
# This function identifies columns (called beads due to old luminex column name nomenclature), given a user-defined percentage of missing values.

# Input ####
# mat: data matrix
# percent.cutoff: numeric; percentage of values with NA values (e.g. for columns with 30% NA values, place 30)

# Output ####
# data frame; listing the column names flagged for removal, and their position in the matrix

beads.rm <- function(mat, percent.cutoff){
  beads.na <- c()
  beads <- colnames(mat)
  for(i in 1:ncol(mat)){
    beads.na[i] <- 100*length(which(is.na(mat[,i])))/length(mat[,i])
  }
  names(beads.na) <- beads
  beads.rm <- names(which(beads.na > percent.cutoff))
  beads.rm.pos <- which(colnames(mat) %in% beads.rm)
  df <- data.frame(beads.rm, beads.rm.pos)
  df
}  # creates a list of beads that were removed

#### Remove columns ####

# This function removes columns (called beads due to old luminex column name nomenclature), given a user-defined percentage of missing values.

# Input ####
# a: data matrix
# percent.cutoff: numeric; percentage of values with NA values (e.g. for columns with 30% NA values, place 30)

# Output ####
# matrix:  columns flagged for removal are removed

remove.beads <- function(a, percent.cutoff){
  beads.na <- c()
  a$remove <- NA
  beads <- colnames(a)
  for(i in 1:ncol(a)){
    beads.na[i] <- 100*length(which(is.na(a[,i])))/length(a[,i])
  }
  names(beads.na) <- beads
  beads.rm <- names(which(beads.na > percent.cutoff))
  beads.rm.pos <- which(colnames(a) %in% beads.rm)
  a <- a[,-c(beads.rm.pos)]
  a
} # removes beads

# Remove Low Counts ####

# This function removes columns (features) that have low representation given a specific value.  This is useful if you want to remove features that are not only present in a certain threshold of samples, but present by at least a certain value

# Input ####
# a: data matrix
# abundance.cutoff: the "zero value" to filter by
# percent.cutoff: the percentage of samples above which the feature should be present in

# example, to select fectures present as a value of 3 in at least 5% of samples, remove.low.counts(a, abundance.cutoff = 3, percent.cutoff = 5)

# Output ####
# list; df = data frame with low abundance features removed, otus.keep = character vector of all features retained

remove.low.counts <- function(a, abundance.cutoff, percent.cutoff){

  beads.low <- c()
  a$remove <- NA
  beads <- colnames(a)
  for(i in 1:ncol(a)){
    beads.low[i] = 100*length(which(a[,i] >= abundance.cutoff )) / length(a[,i])
  }
  names(beads.low) <- beads
  beads.rm <- names(which(beads.low <= percent.cutoff))
  beads.rm <- c(beads.rm, "remove")
  beads.rm.pos <- which(colnames(a) %in% beads.rm)
  a <- a[,-c(beads.rm.pos)]
  data <- list(a, ncol(a))
  names(data) <- c("df", "otus.keep")
  return(data)
}


#### List rows to remove ####
# This function identifies rows (called sub as these are commonly subjects), given a user-defined percentage of missing values.

# Input ####
# mat: data matrix
# percent.cutoff: numeric; percentage of values with NA values (e.g. for columns with 30% NA values, place 30)

# Output ####
# data frame; listing the row names flagged for removal, and their position in the matrix

which.sub.rm <- function(a, percent.cutoff){
  subj <- rownames(a)
  sub.na <- c()
  
  for(i in 1:nrow(a)){
    sub.na[i] <- 100*length(which(is.na(a[i,])))/length(a[i,])
  }
  names(sub.na) <- subj
  subj.rm <- names(which(sub.na > percent.cutoff))
  subj.rm.pos <- which(rownames(a) %in% subj.rm)
  df <- data.frame(subj.rm, subj.rm.pos)
  df
}

#### Remove rows ####

# This function removes rows (called sub as these are commonly subjects), given a user-defined percentage of missing values.

# Input ####
# a: data matrix
# percent.cutoff: numeric; percentage of values with NA values (e.g. for columns with 30% NA values, place 30)

# Output ####
# matrix:  rows flagged for removal are removed

subj.rm <- function(a, percent.cutoff){
  sub.na <- c()
  a[nrow(a) + 1,] <- NA
  subj <- rownames(a)
  for(i in 1:nrow(a)){
    sub.na[i] <- 100*length(which(is.na(a[i,])))/length(a[i,])
  }
  names(sub.na) <- subj
  subj.rm <- names(which(sub.na > 30))
  subj.rm.pos <- which(rownames(a) %in% subj.rm)
  a <- a[-c(subj.rm.pos),]
  a
}

#### Replace NAs ####
# this function replaces missing values by the mean value for that column

# Input: data frame
# Output: data matrix, with missing values replaced by column mean

na.replace <- function(df){
  f1 <- function(vec) { 
    m <- mean(vec, na.rm = TRUE) 
    vec[is.na(vec)] <- m 
    return(vec) 
  } 
  df.narm <- apply(df, 2, f1)
  colnames(df.narm) <- colnames(df)
  rownames(df.narm) <- rownames(df)
  df.narm <- as.matrix(df.narm)
} # replaces NA with column mean