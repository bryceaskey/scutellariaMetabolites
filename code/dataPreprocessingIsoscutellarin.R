library(tidyverse)

rawData <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/raw_data/20190813_fresh.csv")
wrightii <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/raw_data/20210119_freshWrightii.csv")[22:30, ]
rawData <- rbind(rawData, wrightii, make.row.names=FALSE)

# Define functions for interpreting injection names
abbrevNames <- data.frame(
  abbrev=c("HV", "AC", "AS", "BT", "TT", "HF", "BL", "LD", "RMSEQ", "R071119", "R_MS", "R_SC", "WRI"),
  fullName=c("havanensis", "arenicola", "altissima", "barbata", "tournefortii", "hastifolia", "baicalensis", "leonardii", "racemosa_RNAseq", "racemosa_071119", "racemosa_MS", "racemosa_SC", "wrightii")
)

getSampleName <- function(injectionName){
  sampleAbbrev <- strsplit(injectionName, "-")[[1]][1]
  sampleName <- abbrevNames[grep(sampleAbbrev, abbrevNames$abbrev), 2]
  return(sampleName)
}

getSampleRep <- function(injectionName){
  sampleRep <- strsplit(injectionName, "-")[[1]][2]
  return(sampleRep)
}

abbrevOrgans <- data.frame(
  abbrev=c("L", "S", "R"),
  fullOrgan=c("leaves", "stems", "roots")
)

getSampleOrgan <- function(injectionName){
  organAbbrev <- strsplit(injectionName, "-")[[1]][3]
  sampleOrgan <- abbrevOrgans[grep(organAbbrev, abbrevOrgans$abbrev), 2]
  return(sampleOrgan)
}

rawDataSubset <- data.frame(injectionName=rawData$injectionName,
                            isoscutellarin=rawData$isoscutellarin)

# Create dataframe with preprocessed data
dataRows <- c(23:112, 125:133)
allData <- data.frame(
  species=sapply(rawData$injectionName[dataRows], getSampleName),
  replicate=sapply(rawData$injectionName[dataRows], getSampleRep),
  organ=sapply(rawData$injectionName[dataRows], getSampleOrgan),
  isoscutellarin=rawDataSubset$isoscutellarin[dataRows]
)

# Calculate mean and standard error
allData <- pivot_longer(allData, cols=c(4), names_to="metabolite", values_to="peakArea")
allData <- allData %>%
  group_by(species, organ, metabolite) %>%
  summarise(mean_peakArea=mean(peakArea), stError_peakArea=sd(peakArea)/sqrt(3))
colnames(allData)[4] <- "peakArea"

write.csv(allData, file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20190813_isoscutellarin.csv", row.names=FALSE)
