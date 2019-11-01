# Injection names should be standardized to match the following format:
# Variety Abbreviation - Replicate Number - Plant Organ
# e.g. BL-1-R = baicalensis - 1st replicate - Root

library(ggplot2)
library(plyr)
library(tibble)

# Read metabolite data from .csv file -------------------------------------------------------------
rawData <- read.csv(file="C:/Users/bca08_000/Documents/scutellaria/data/metaboliteData.csv", header=TRUE)
rawData[, 1] <- as.character(rawData[, 1])

# Define functions for interpreting injection names -----------------------------------------------
abbrevNames <- data.frame(
  abbrev=c("HV", "AC", "AS", "BT", "TT", "HF", "BL", "LD", "RMSEQ", "R071119", "R_MS", "R_SC"),
  fullName=c("Havenesis", "Arenicola", "Altissima", "Barbata", "Tourmetii", "Hastafolia", "Baicalensis", "Lateriflora", "RNA Seq", "Racemosa 071119", "Racemosa MS", "Racemosa SC")
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
  fullOrgan=c("Leaves", "Shoots", "Roots")
)

getSampleOrgan <- function(injectionName){
  organAbbrev <- strsplit(injectionName, "-")[[1]][3]
  sampleOrgan <- abbrevOrgans[grep(organAbbrev, abbrevOrgans$abbrev), 2]
  return(sampleOrgan)
}

# Define function to convert peak area data into ppm using calibration samples --------------------
ppmConversion <- function(peakArea, metaboliteName, rawData=rawData){
  if(sum(grepl(metaboliteName, colnames(rawData)[2:9]))==1){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][5:12])
  }else{
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][14:21])
  }
  convFactor <- lm(formula=ppm~area+0, data=calibrations)[[1]][[1]]
  ppm <- peakArea*convFactor
  return(ppm)
}

# Create data frame with processed data -----------------------------------------------------------
allData <- data.frame(
  variety=sapply(rawData$injectionName[23:112], getSampleName),
  replicate=sapply(rawData$injectionName[23:112], getSampleRep),
  organ=sapply(rawData$injectionName[23:112], getSampleOrgan),
  apigenin=sapply(rawData$apigenin[23:112], ppmConversion, metaboliteName="apigenin", rawData=rawData),
  apigeninG=sapply(rawData$apigeninG[23:112], ppmConversion, metaboliteName="apigeninG", rawData=rawData),
  scutellarein=sapply(rawData$scutellarein[23:112], ppmConversion, metaboliteName="scutellarein", rawData=rawData),
  scutellarin=sapply(rawData$scutellarin[23:112], ppmConversion, metaboliteName="scutellarin", rawData=rawData),
  hispidulin=sapply(rawData$hispidulin[23:112], ppmConversion, metaboliteName="hispidulin", rawData=rawData),
  hispiduloG=sapply(rawData$hispiduloG[23:112], ppmConversion, metaboliteName="hispiduloG", rawData=rawData),
  chrysin=sapply(rawData$chrysin[23:112], ppmConversion, metaboliteName="chrysin", rawData=rawData),
  chrysinG=sapply(rawData$chrysinG[23:112], ppmConversion, metaboliteName="chrysinG", rawData=rawData),
  wogonin=sapply(rawData$wogonin[23:112], ppmConversion, metaboliteName="wogonin", rawData=rawData),
  wogonoside=sapply(rawData$wogonoside[23:112], ppmConversion, metaboliteName="wogonoside", rawData=rawData),
  baicalein=sapply(rawData$baicalein[23:112], ppmConversion, metaboliteName="baicalein", rawData=rawData),
  baicalin=sapply(rawData$baicalin[23:112], ppmConversion, metaboliteName="baicalin", rawData=rawData),
  oroxylinA=sapply(rawData$oroxylinA[23:112], ppmConversion, metaboliteName="oroxylinA", rawData=rawData),
  oroxyloside=sapply(rawData$oroxyloside[23:112], ppmConversion, metaboliteName="oroxyloside", rawData=rawData),
  acetoside=sapply(rawData$acetoside[23:112], ppmConversion, metaboliteName="acetoside", rawData=rawData)
)

# Calculate means and standard deviations for each metabolite in each variety ---------------------
meanData <- ddply(allData, c("Variety", "Organ"), summarise,
  acetosideMean=mean(...),
  acetosideStError=sd(...)/sqrt(length(...))
)