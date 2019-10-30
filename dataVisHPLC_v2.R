library(ggplot2)
library(plyr)
library(tibble)

# Read metabolite data csv file
rawData <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/metaboliteData.csv", header=TRUE)
metaboliteData <- rawData[23:112, ]
metaboliteData[, 1] <- as.character(metaboliteData[, 1])

# Define functions for interpreting injection names
abbrevNames <- data.frame(
  abbrev=c("HV", "AC", "AS", "BT", "TT", "HF", "BL", "LD", "RMSEQ", "R(MS)", "R(SC)"),
  fullName=c("Havenesis", "Arenicola", "Altissima", "Barbata", "Tourmetii", "Hastafolia", "Baicalensis", "Lateriflora", "RNA Seq", "Racemosa MS", "Racemosa SC")
)

getSampleName <- function(injectionName){
  sampleAbbrev <- strsplit(injectionName, "-")[[1]][1]
  sampleName <- abbrevNames[grep(sampleAbbrev, abbrevNames$abbrev), 2]
  return(sampleName)
}

abbrevOrgans <- data.frame(
  abbrev=c("L", "S", "R"),
  fullOrgan=c("Leaves", "Shoots", "Roots")
)

getSampleOrgan <- function(injectionName){
  organAbbrev <- strsplit(injectionName, "-")[[1]][3]
  sampleOrgan <- abbrevOrgans[grep(organAbbrev, abbrevOrgans$abbrev), 2]
  return(sampleOrgan)
}



#for(row in 1:nrow(metaboliteData)){
#  sampleName <- append(sampleName, strsplit(metaboliteData$injectionName[row], "-")[[1]][1])
#  plantOrgan <- append(plantOrgan, strsplit(metaboliteData$injectionName[row], "-")[[1]][3])
  #metReplicates <- rbind(metaboliteData[row, 2:16], metaboliteData[row+1, 2:16], metaboliteData[row+2, 2:16])
  #avgMetData <- rbind(avgMetData, colMeans(metReplicates))
  #metabolites <- colMeans(subset(metaboliteData, grepl(metaboliteData$injectionName, sampleName[length(sampleName)]), select=2:16))
#}