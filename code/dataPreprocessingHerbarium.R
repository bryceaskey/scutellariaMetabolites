# Data preprocessing for herbarium HPLC data
# Finds species name from injection #, uses calibration standards to calculate ppm from peak area, and shapes data into tidy format
library(tidyr)

rawHerbariumData1_30 <- read.csv(file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/raw_data/20191217_dryHerbarium1-30.csv", header=TRUE)
rawHerbariumData31_78 <- read.csv(file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/raw_data/20200714_dryHerbarium31-78.csv", header=TRUE)
herbariumList <- read.csv(file="C:/Users/Bryce/Research/scutellariaMetabolites/data/herbarium/herbariumList.csv", header=TRUE)[, 1:4]

getSampleName <- function(injectionName, herbariumList){
  sampleName <- herbariumList[injectionName==herbariumList$Label, 3]
  if(length(sampleName)==0){
    sampleName <- injectionName
    print(paste("Warning, injection name:", injectionName, "not recognized"))
  }
  return(as.character(sampleName))
}

getSampleOrgan <- function(injectionName, herbariumList){
  sampleOrgan <- herbariumList[injectionName==herbariumList$Label, 4]
  if(length(sampleOrgan)==0){
    sampleOrgan <- NA
  }
  return(as.character(sampleOrgan))
}

# Function for samples 1-30
ppmConversion <- function(peakArea, metaboliteName, rawData=rawData){
  mix1Metabolites <- c("acteoside", "baicalin", "oroxyloside", "wogonoside", "apigenin", "baicalein", "wogonin", "oroxylinA")
  mix2Metabolites <- c("hispidulinG", "apigeninG", "isoscutellarin", "scutellarein", "chrysinG", "hispidulin", "chrysin")
  mix1Rows <- c(11:18)
  mix2Rows <- c(20:27)
  scutellarinRows <- c(75:82)
  
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

dataRows <- c(30:73)
herbariumData1_30 <- data.frame(
  species=sapply(rawHerbariumData1_30$injection[dataRows], getSampleName, herbariumList=herbariumList),
  organ=sapply(rawHerbariumData1_30$injection[dataRows], getSampleOrgan, herbariumList=herbariumList),
  acteoside=sapply(rawHerbariumData1_30$acteoside[dataRows], ppmConversion, metaboliteName="acteoside", rawData=rawHerbariumData1_30),
  baicalin=sapply(rawHerbariumData1_30$baicalin[dataRows], ppmConversion, metaboliteName="baicalin", rawData=rawHerbariumData1_30),
  oroxyloside=sapply(rawHerbariumData1_30$oroxyloside[dataRows], ppmConversion, metaboliteName="oroxyloside", rawData=rawHerbariumData1_30),
  wogonoside=sapply(rawHerbariumData1_30$wogonoside[dataRows], ppmConversion, metaboliteName="wogonoside", rawData=rawHerbariumData1_30),
  apigenin=sapply(rawHerbariumData1_30$apigenin[dataRows], ppmConversion, metaboliteName="apigenin", rawData=rawHerbariumData1_30),
  baicalein=sapply(rawHerbariumData1_30$baicalein[dataRows], ppmConversion, metaboliteName="baicalein", rawData=rawHerbariumData1_30),
  wogonin=sapply(rawHerbariumData1_30$wogonin[dataRows], ppmConversion, metaboliteName="wogonin", rawData=rawHerbariumData1_30),
  oroxylinA=sapply(rawHerbariumData1_30$oroxylinA[dataRows], ppmConversion, metaboliteName="oroxylinA", rawData=rawHerbariumData1_30),
  hispidulinG=sapply(rawHerbariumData1_30$hispidulinG[dataRows], ppmConversion, metaboliteName="hispidulinG", rawData=rawHerbariumData1_30),
  apigeninG=sapply(rawHerbariumData1_30$apigeninG[dataRows], ppmConversion, metaboliteName="apigeninG", rawData=rawHerbariumData1_30),
  isoscutellarin=sapply(rawHerbariumData1_30$isoscutellarin[dataRows], ppmConversion, metaboliteName="isoscutellarin", rawData=rawHerbariumData1_30),
  scutellarein=sapply(rawHerbariumData1_30$scutellarein[dataRows], ppmConversion, metaboliteName="scutellarein", rawData=rawHerbariumData1_30),
  chrysinG=sapply(rawHerbariumData1_30$chrysinG[dataRows], ppmConversion, metaboliteName="chrysinG", rawData=rawHerbariumData1_30),
  hispidulin=sapply(rawHerbariumData1_30$hispidulin[dataRows], ppmConversion, metaboliteName="hispidulin", rawData=rawHerbariumData1_30),
  chrysin=sapply(rawHerbariumData1_30$chrysin[dataRows], ppmConversion, metaboliteName="chrysin", rawData=rawHerbariumData1_30),
  scutellarin=sapply(rawHerbariumData1_30$scutellarin[dataRows], ppmConversion, metaboliteName="scutellarin", rawData=rawHerbariumData1_30)
)

# Function for samples 31-78
ppmConversion <- function(peakArea, metaboliteName, rawData=rawHerbariumData31_78){
  mix1Metabolites <- c("acteoside", "baicalin", "oroxyloside", "wogonoside", "apigenin", "baicalein", "wogonin", "oroxylinA")
  mix2Metabolites <- c("hispidulinG", "apigeninG", "isoscutellarin", "scutellarein", "chrysinG", "hispidulin", "chrysin")
  mix1Rows <- c(5:11)
  mix2Rows <- c(14:20)
  scutellarinRows <- c(72:79)
  
  if(metaboliteName %in% mix1Metabolites){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50), area=rawData[[metaboliteName]][mix1Rows])
  }else if(metaboliteName %in% mix2Metabolites){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50), area=rawData[[metaboliteName]][mix2Rows])
  }else if(metaboliteName == "scutellarin"){
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][scutellarinRows])
  }else{
    stop(paste("Error: metabolite name not recognized:", metaboliteName))
  }
  
  convFactor <- lm(formula=ppm~area+0, data=calibrations)[[1]][[1]]
  ppm <- peakArea*convFactor
  return(ppm)
}



dataRows <- c(23:70)
herbariumData31_78 <- data.frame(
  species=sapply(rawHerbariumData31_78$injection[dataRows], getSampleName, herbariumList=herbariumList),
  organ=sapply(rawHerbariumData31_78$injection[dataRows], getSampleOrgan, herbariumList=herbariumList),
  acteoside=sapply(rawHerbariumData31_78$acteoside[dataRows], ppmConversion, metaboliteName="acteoside", rawData=rawHerbariumData31_78),
  baicalin=sapply(rawHerbariumData31_78$baicalin[dataRows], ppmConversion, metaboliteName="baicalin", rawData=rawHerbariumData31_78),
  oroxyloside=sapply(rawHerbariumData31_78$oroxyloside[dataRows], ppmConversion, metaboliteName="oroxyloside", rawData=rawHerbariumData31_78),
  wogonoside=sapply(rawHerbariumData31_78$wogonoside[dataRows], ppmConversion, metaboliteName="wogonoside", rawData=rawHerbariumData31_78),
  apigenin=sapply(rawHerbariumData31_78$apigenin[dataRows], ppmConversion, metaboliteName="apigenin", rawData=rawHerbariumData31_78),
  baicalein=sapply(rawHerbariumData31_78$baicalein[dataRows], ppmConversion, metaboliteName="baicalein", rawData=rawHerbariumData31_78),
  wogonin=sapply(rawHerbariumData31_78$wogonin[dataRows], ppmConversion, metaboliteName="wogonin", rawData=rawHerbariumData31_78),
  oroxylinA=sapply(rawHerbariumData31_78$oroxylinA[dataRows], ppmConversion, metaboliteName="oroxylinA", rawData=rawHerbariumData31_78),
  hispidulinG=sapply(rawHerbariumData31_78$hispidulinG[dataRows], ppmConversion, metaboliteName="hispidulinG", rawData=rawHerbariumData31_78),
  apigeninG=sapply(rawHerbariumData31_78$apigeninG[dataRows], ppmConversion, metaboliteName="apigeninG", rawData=rawHerbariumData31_78),
  isoscutellarin=sapply(rawHerbariumData31_78$isoscutellarin[dataRows], ppmConversion, metaboliteName="isoscutellarin", rawData=rawHerbariumData31_78),
  scutellarein=sapply(rawHerbariumData31_78$scutellarein[dataRows], ppmConversion, metaboliteName="scutellarein", rawData=rawHerbariumData31_78),
  chrysinG=sapply(rawHerbariumData31_78$chrysinG[dataRows], ppmConversion, metaboliteName="chrysinG", rawData=rawHerbariumData31_78),
  hispidulin=sapply(rawHerbariumData31_78$hispidulin[dataRows], ppmConversion, metaboliteName="hispidulin", rawData=rawHerbariumData31_78),
  chrysin=sapply(rawHerbariumData31_78$chrysin[dataRows], ppmConversion, metaboliteName="chrysin", rawData=rawHerbariumData31_78),
  scutellarin=sapply(rawHerbariumData31_78$scutellarin[dataRows], ppmConversion, metaboliteName="scutellarin", rawData=rawHerbariumData31_78)
)

speciesToRemove <- paste(c("15mix", "wash"), collapse = '|')
herbariumData1_30 <- filter(herbariumData1_30, !grepl(speciesToRemove, species))
herbariumData31_78 <- filter(herbariumData31_78, !grepl(speciesToRemove, species))

herbariumData1_30 <- pivot_longer(herbariumData1_30, cols=c(3:18), names_to="metabolite", values_to="concentration_ppm")
herbariumData1_30$stError_ppm <- NA

herbariumData31_78 <- pivot_longer(herbariumData31_78, cols=c(3:18), names_to="metabolite", values_to="concentration_ppm")
herbariumData31_78$stError_ppm <- NA

write.csv(herbariumData1_30, file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200214_herbarium1_30.csv", row.names=FALSE)
write.csv(herbariumData31_78, file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200812_herbarium31_78.csv", row.names=FALSE)