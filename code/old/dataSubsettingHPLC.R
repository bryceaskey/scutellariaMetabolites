# Function to subset data frame output by dataExtraction function into smaller data frames.
# Useful if data from multiple experiments is mixed together, or if wash data is mixed in with experimental data.

dataSubsetting <- function(inputDataFrame, mode="remove", sample, rowInds=NULL)
removeInds <- c(grep("JA", rownames(areaData)), grep("Ag