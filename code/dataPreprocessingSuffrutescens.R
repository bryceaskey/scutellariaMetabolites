library(tidyverse)

rawData <- read.csv(file="C:/Users/bca08_000/Documents/scutellariaMetabolites/data/20201119_suffrutescens.csv", header=TRUE)
colnames(rawData)[1] <- "injectionName"

# Define function to interpret injection names ----
getSampleName <- function(injectionName){
  sampleName <- "suffrutescens"
  return(as.character(sampleName))
}

# Define function to interpret injection organ ----
getSampleOrgan <- function(injectionName){
  organAbbrev <- strsplit(strsplit(injectionName, split="_")[[1]][2], split="-")[[1]][1]
  if(organAbbrev=="L"){
    sampleOrgan <- "leaves"
  }else if(organAbbrev=="S"){
    sampleOrgan <- "shoots"
  }else if(organAbbrev=="R"){
    sampleOrgan <- "roots"
  }else{
    sampleOrgan <- NA
  }
  return(as.character(sampleOrgan))
}

# Define function to convert peak area data into ppm using calibration samples ----
ppmConversion <- function(peakArea, flavonoidName, mix1Rows, mix2Rows, rawData){
  mix1Flavonoids <- c("acetoside", "baicalin", "oroxyloside", "wogonoside", "apigenin", "baicalein", "wogonin", "oroxylinA")
  mix2Flavonoids <- c("hispidulinG", "apigeninG", "scutellarin", "scutellarein", "chrysinG", "hispidulin", "chrysin")
  if(length(mix1Rows)==8){
    ppmList <- c(0.1, 0.5, 1, 5, 10, 25, 50, 100)
  }else{
    ppmList <- c(0.1, 0.5, 1, 5, 10, 25, 50)
  }
  
  if(flavonoidName %in% mix1Flavonoids){
    calibrations <- data.frame(ppm=ppmList, area=rawData[[flavonoidName]][mix1Rows])
    convFactor <- lm(formula=ppm~area+0, data=calibrations)[[1]][[1]]
    ppm <- peakArea*convFactor
    return(ppm)
  }else if(flavonoidName %in%mix2Flavonoids){
    calibrations <- data.frame(ppm=ppmList, area=rawData[[flavonoidName]][mix2Rows])
    convFactor <- lm(formula=ppm~area+0, data=calibrations)[[1]][[1]]
    ppm <- peakArea*convFactor
    return(ppm)
  }else{
    print("Flavonoid name not recognized")
    return(NA)
  }
}

mix1Rows <- c(4:11)
mix2Rows <- c(13:20)
dataRows <- c(22:30)
sufData <- data.frame(
  species=sapply(rawData$injection[dataRows], getSampleName),
  organ=sapply(rawData$injection[dataRows], getSampleOrgan),
  apigenin=sapply(rawData$apigenin[dataRows], ppmConversion, flavonoidName="apigenin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  apigeninG=sapply(rawData$apigeninG[dataRows], ppmConversion, flavonoidName="apigeninG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  scutellarein=sapply(rawData$scutellarein[dataRows], ppmConversion, flavonoidName="scutellarein", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  scutellarin=sapply(rawData$scutellarin[dataRows], ppmConversion, flavonoidName="scutellarin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  hispidulin=sapply(rawData$hispidulin[dataRows], ppmConversion, flavonoidName="hispidulin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  hispidulinG=sapply(rawData$hispidulinG[dataRows], ppmConversion, flavonoidName="hispidulinG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  chrysin=sapply(rawData$chrysin[dataRows], ppmConversion, flavonoidName="chrysin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  chrysinG=sapply(rawData$chrysinG[dataRows], ppmConversion, flavonoidName="chrysinG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  wogonin=sapply(rawData$wogonin[dataRows], ppmConversion, flavonoidName="wogonin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  wogonoside=sapply(rawData$wogonoside[dataRows], ppmConversion, flavonoidName="wogonoside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  baicalein=sapply(rawData$baicalein[dataRows], ppmConversion, flavonoidName="baicalein", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  baicalin=sapply(rawData$baicalin[dataRows], ppmConversion, flavonoidName="baicalin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  oroxylinA=sapply(rawData$oroxylinA[dataRows], ppmConversion, flavonoidName="oroxylinA", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  oroxyloside=sapply(rawData$oroxyloside[dataRows], ppmConversion, flavonoidName="oroxyloside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData),
  acetoside=sapply(rawData$acetoside[dataRows], ppmConversion, flavonoidName="acteoside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawData)
)

sufData <- pivot_longer(sufData, cols=c(3:17), names_to="metabolite", values_to="concentration_ppm")
sufData <- sufData %>%
  group_by(species, organ, metabolite) %>%
  summarise(concentration_ppm_mean=mean(concentration_ppm), stError_ppm=sd(concentration_ppm)/sqrt(3))
colnames(sufData)[4] <- "concentration_ppm"

write.csv(sufData, file="C:/Users/bca08_000/Documents/scutellariaMetabolites/data/preprocessed/20201119_suffrutescens.csv")
