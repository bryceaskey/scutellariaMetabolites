library(tidyverse)

rawData <- read.csv(file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/raw_data/20200924_freshWrightii.csv", header=TRUE)

# Define function to interpret injection names
getSampleName <- function(injectionName){
  sampleName <- "wrightii"
  return(as.character(sampleName))
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

# Define function to convert peak area data into ppm using calibration samples
ppmConversion <- function(peakArea, metaboliteName, rawData=rawData){
  mix1Metabolites <- c("acteoside", "baicalin", "oroxyloside", "wogonoside", "apigenin", "baicalein", "wogonin", "oroxylinA")
  mix2Metabolites <- c("hispidulinG", "apigeninG", "isoscutellarin", "scutellarein", "chrysinG", "hispidulin", "chrysin")
  mix1Rows <- 5:12
  mix2Rows <- 14:21
  scutellarinRows <- 35:42
  
  if(metaboliteName %in% mix1Metabolites){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][mix1Rows])
  }else if(metaboliteName %in% mix2Metabolites){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][mix2Rows])
  }else if(metaboliteName == "scutellarin"){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][scutellarinRows])
  }else{
    stop(paste("Error: metabolite name not recognized:", metaboliteName))
  }
  
  convFactor <- lm(formula=ppm~area+0, data=calibrations)[[1]][[1]]
  ppm <- peakArea*convFactor
  return(ppm)
}

# Create dataframe with preprocessed data
dataRows <- c(23:31)
allData <- data.frame(
  species=sapply(rawData$injectionName[dataRows], getSampleName),
  replicate=sapply(rawData$injectionName[dataRows], getSampleRep),
  organ=sapply(rawData$injectionName[dataRows], getSampleOrgan),
  acteoside=sapply(rawData$acteoside[dataRows], ppmConversion, metaboliteName="acteoside", rawData=rawData),
  baicalin=sapply(rawData$baicalin[dataRows], ppmConversion, metaboliteName="baicalin", rawData=rawData),
  oroxyloside=sapply(rawData$oroxyloside[dataRows], ppmConversion, metaboliteName="oroxyloside", rawData=rawData),
  wogonoside=sapply(rawData$wogonoside[dataRows], ppmConversion, metaboliteName="wogonoside", rawData=rawData),
  apigenin=sapply(rawData$apigenin[dataRows], ppmConversion, metaboliteName="apigenin", rawData=rawData),
  baicalein=sapply(rawData$baicalein[dataRows], ppmConversion, metaboliteName="baicalein", rawData=rawData),
  wogonin=sapply(rawData$wogonin[dataRows], ppmConversion, metaboliteName="wogonin", rawData=rawData),
  oroxylinA=sapply(rawData$oroxylinA[dataRows], ppmConversion, metaboliteName="oroxylinA", rawData=rawData),
  hispidulinG=sapply(rawData$hispidulinG[dataRows], ppmConversion, metaboliteName="hispidulinG", rawData=rawData),
  apigeninG=sapply(rawData$apigeninG[dataRows], ppmConversion, metaboliteName="apigeninG", rawData=rawData),
  isoscutellarin=sapply(rawData$isoscutellarin[dataRows], ppmConversion, metaboliteName="isoscutellarin", rawData=rawData),
  scutellarein=sapply(rawData$scutellarein[dataRows], ppmConversion, metaboliteName="scutellarein", rawData=rawData),
  chrysinG=sapply(rawData$chrysinG[dataRows], ppmConversion, metaboliteName="chrysinG", rawData=rawData),
  hispidulin=sapply(rawData$hispidulin[dataRows], ppmConversion, metaboliteName="hispidulin", rawData=rawData),
  chrysin=sapply(rawData$chrysin[dataRows], ppmConversion, metaboliteName="chrysin", rawData=rawData),
  scutellarin=sapply(rawData$scutellarin[dataRows], ppmConversion, metaboliteName="scutellarin", rawData=rawData)
)

allData <- pivot_longer(allData, cols=c(4:19), names_to="metabolite", values_to="concentration_ppm")
allData <- allData %>%
  group_by(species, organ, metabolite) %>%
  summarise(concentration_ppm_mean=mean(concentration_ppm), stError_ppm=sd(concentration_ppm)/sqrt(3))
colnames(allData)[4] <- "concentration_ppm"

write.csv(allData, file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201007_wrightii.csv", row.names=FALSE)
