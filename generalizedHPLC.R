# Script to extract peak height and area data from HPLC excel files exported by Chromeleon software.
# Note - to avoid confusion regarding peak identity, only data from named peaks will be saved and analyzed.

# install any necessary packages if they aren't already installed
packageList <- c("readxl", "ggplot2")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0){
  install.packages(newPackages)
}

# load packages into environment
library(readxl)
library(ggplot2)

# set folder containing all HPLC excel files as the working directory
setwd("C:/Users/bca08_000/Documents/scutellaria/Test2 - Processing 2")

# initialize empty list to store dataframes with all integration data
allData <- list()

# loop through all files specified by setwd()
for(i in 1:length(dir())){
  # read excel file into R as a dataframe - only sheet named "Integration"
  fullSheet <- read_excel(dir()[i], sheet="Integration", .name_repair="minimal")
  # extract section of fullSheet dataframe containing integration results
  startRow <- which(fullSheet[ , 1] == "Integration Results") + 4
  endRow <- nrow(fullSheet) - 1
  intResults <- fullSheet[startRow:endRow, ]
  # assign column names
  colnames(intResults) <- paste(fullSheet[startRow - 3, ])
  # correct formatting/syntax issuses, and remove unnamed peaks
  intResults[intResults=="n.a."] <- NA
  intResults <- intResults[is.na(intResults$`Peak Name`)==FALSE, ]
  intResults[, 3:ncol(intResults)] <- lapply(3:ncol(intResults), function(x) as.numeric(intResults[[x]]))
  # save dataframe into list, and name
  allData[[i]] <- intResults
  names(allData)[i] <- dir()[i]
}

# TODO: add optional method to correct for dilutions - standardize format for dilution files
# can dilution factor be set in Chromeleon software?

# TODO: use subset function to identify and separate experiemental data from standard/wash data in allData
# or have user separate files manually?

# create new data frames to store peak height and area data for all samples
areaData <- data.frame(matrix(ncol=length(intResults$`Peak Name`), nrow=length(allData)))
heightData <- data.frame(matrix(ncol=length(intResults$`Peak Name`), nrow=length(allData)))

# name columns of height and area dataframes to match peak names
colnames(areaData) <- intResults$`Peak Name`
colnames(heightData) <- intResults$`Peak Name`

# loop through each dataframe in allData
for(i in 1:length(allData)){
  # save area and height data for each named peak
  areaData[i, ] <- allData[[i]]$Area
  heightData[i, ] <- allData[[i]]$Height 
  # name row to match the name of the sample
  rownames(areaData)[i] <- names(allData[i])
  rownames(heightData)[i] <- names(allData[i])
}



