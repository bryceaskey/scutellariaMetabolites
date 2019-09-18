# install necessary packages if they aren't already installed
packageList <- c("readxl", "ggplot2")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0){
  install.packages(newPackages)
}

# load packages into environment
library(readxl)
library(ggplot2)

# read HPLC data from excel files, and store as a list of dataframes
setwd("C:/Users/bca08_000/Documents/scutellaria/Test2 - Processing 2") # folder containing all HPLC excel files
allData <- list()
for(i in 1:length(dir())){
  fullSheet <- read_excel(dir()[i], sheet="Integration", .name_repair="minimal")
  startRow <- which(fullSheet[ , 1] == "Integration Results") + 4
  endRow <- nrow(fullSheet) - 1
  intResults <- fullSheet[startRow:endRow, ]
  # don't include n.a. peaks
  colnames(intResults) <- paste(fullSheet[startRow - 3, ])
  intResults[intResults=="n.a."] <- NA
  allData[[i]] <- intResults
}
names(allData) <- dir()

# need method to correct for dilutions

# 

