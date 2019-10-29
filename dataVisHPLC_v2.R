library(ggplot2)
library(plyr)

rawData <- read.csv(file="C:/Users/bca08_000/Documents/scutellaria/metaboliteData.csv", header=TRUE)
metaboliteData <- rawData[23:112, ]

abbreviatedNames <- data.frame(
  abbreviation=c("HV", "AC", "AS", "BT", "TT", "HF", "BL", "LD"),
  fullName=c("Havenesis", "Arenicola", "Altissima", "Barbata", "Tourmetii", "Hastafolia", "Baicalensis", "Lateriflora")
)

getSampleName <- function(injectionName){
  sampleAbbrev <- strsplit(injectionName, "-")[[1]][1]
  
}

#for(row in 1:nrow(metaboliteData)){
#  sampleName <- append(sampleName, strsplit(metaboliteData$injectionName[row], "-")[[1]][1])
#  plantOrgan <- append(plantOrgan, strsplit(metaboliteData$injectionName[row], "-")[[1]][3])
  #metReplicates <- rbind(metaboliteData[row, 2:16], metaboliteData[row+1, 2:16], metaboliteData[row+2, 2:16])
  #avgMetData <- rbind(avgMetData, colMeans(metReplicates))
  #metabolites <- colMeans(subset(metaboliteData, grepl(metaboliteData$injectionName, sampleName[length(sampleName)]), select=2:16))
#}