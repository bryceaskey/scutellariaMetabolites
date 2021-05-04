# Data preprocessing for herbarium HPLC data
# Finds species name from injection #, uses calibration standards to calculate ppm from peak area, and shapes data into tidy format
library(tidyr)

rawHerbariumData1_30 <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/metaboliteDataHerbarium1-30.csv", header=TRUE)
rawHerbariumData1_30 <- rawHerbariumData1_30[1:78,]
rawHerbariumData31_78 <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/metaboliteDataHerbarium31-78.csv", header=TRUE)
herbariumList <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/herbariumList.csv", header=TRUE)

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

# Define function to convert peak area data into ppm using calibration samples --------------------
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

mix1Rows <- c(11:18)
mix2Rows <- c(20:27)
dataRows <- c(28:78)
herbariumData1_30 <- data.frame(
  species=sapply(rawHerbariumData1_30$injection[dataRows], getSampleName, herbariumList=herbariumList),
  organ=sapply(rawHerbariumData1_30$injection[dataRows], getSampleOrgan, herbariumList=herbariumList),
  apigenin=sapply(rawHerbariumData1_30$apigenin[dataRows], ppmConversion, flavonoidName="apigenin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  apigeninG=sapply(rawHerbariumData1_30$apigeninG[dataRows], ppmConversion, flavonoidName="apigeninG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  scutellarein=sapply(rawHerbariumData1_30$scutellarein[dataRows], ppmConversion, flavonoidName="scutellarein", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  scutellarin=sapply(rawHerbariumData1_30$scutellarin[dataRows], ppmConversion, flavonoidName="scutellarin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  hispidulin=sapply(rawHerbariumData1_30$hispidulin[dataRows], ppmConversion, flavonoidName="hispidulin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  hispidulinG=sapply(rawHerbariumData1_30$hispidulinG[dataRows], ppmConversion, flavonoidName="hispidulinG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  chrysin=sapply(rawHerbariumData1_30$chrysin[dataRows], ppmConversion, flavonoidName="chrysin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  chrysinG=sapply(rawHerbariumData1_30$chrysinG[dataRows], ppmConversion, flavonoidName="chrysinG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  wogonin=sapply(rawHerbariumData1_30$wogonin[dataRows], ppmConversion, flavonoidName="wogonin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  wogonoside=sapply(rawHerbariumData1_30$wogonoside[dataRows], ppmConversion, flavonoidName="wogonoside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  baicalein=sapply(rawHerbariumData1_30$baicalein[dataRows], ppmConversion, flavonoidName="baicalein", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  baicalin=sapply(rawHerbariumData1_30$baicalin[dataRows], ppmConversion, flavonoidName="baicalin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  oroxylinA=sapply(rawHerbariumData1_30$oroxylinA[dataRows], ppmConversion, flavonoidName="oroxylinA", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  oroxyloside=sapply(rawHerbariumData1_30$oroxyloside[dataRows], ppmConversion, flavonoidName="oroxyloside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30),
  acetoside=sapply(rawHerbariumData1_30$acetoside[dataRows], ppmConversion, flavonoidName="acetoside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData1_30)
)

mix1Rows <- c(5:11)
mix2Rows <- c(14:20)
dataRows <- c(22:75)
herbariumData31_78 <- data.frame(
  species=sapply(rawHerbariumData31_78$injection[dataRows], getSampleName, herbariumList=herbariumList),
  organ=sapply(rawHerbariumData31_78$injection[dataRows], getSampleOrgan, herbariumList=herbariumList),
  apigenin=sapply(rawHerbariumData31_78$apigenin[dataRows], ppmConversion, flavonoidName="apigenin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  apigeninG=sapply(rawHerbariumData31_78$apigeninG[dataRows], ppmConversion, flavonoidName="apigeninG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  scutellarein=sapply(rawHerbariumData31_78$scutellarein[dataRows], ppmConversion, flavonoidName="scutellarein", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  scutellarin=sapply(rawHerbariumData31_78$scutellarin[dataRows], ppmConversion, flavonoidName="scutellarin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  hispidulin=sapply(rawHerbariumData31_78$hispidulin[dataRows], ppmConversion, flavonoidName="hispidulin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  hispidulinG=sapply(rawHerbariumData31_78$hispidulinG[dataRows], ppmConversion, flavonoidName="hispidulinG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  chrysin=sapply(rawHerbariumData31_78$chrysin[dataRows], ppmConversion, flavonoidName="chrysin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  chrysinG=sapply(rawHerbariumData31_78$chrysinG[dataRows], ppmConversion, flavonoidName="chrysinG", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  wogonin=sapply(rawHerbariumData31_78$wogonin[dataRows], ppmConversion, flavonoidName="wogonin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  wogonoside=sapply(rawHerbariumData31_78$wogonoside[dataRows], ppmConversion, flavonoidName="wogonoside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  baicalein=sapply(rawHerbariumData31_78$baicalein[dataRows], ppmConversion, flavonoidName="baicalein", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  baicalin=sapply(rawHerbariumData31_78$baicalin[dataRows], ppmConversion, flavonoidName="baicalin", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  oroxylinA=sapply(rawHerbariumData31_78$oroxylinA[dataRows], ppmConversion, flavonoidName="oroxylinA", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  oroxyloside=sapply(rawHerbariumData31_78$oroxyloside[dataRows], ppmConversion, flavonoidName="oroxyloside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78),
  acetoside=sapply(rawHerbariumData31_78$acetoside[dataRows], ppmConversion, flavonoidName="acetoside", mix1Rows=mix1Rows, mix2Rows=mix2Rows, rawData=rawHerbariumData31_78)
)

speciesToRemove <- paste(c("15mix", "wash"), collapse = '|')
herbariumData1_30 <- filter(herbariumData1_30, !grepl(speciesToRemove, species))
herbariumData31_78 <- filter(herbariumData31_78, !grepl(speciesToRemove, species))

herbariumData1_30 <- pivot_longer(herbariumData1_30, cols=c(3:17), names_to="metabolite", values_to="concentration_ppm")
herbariumData1_30$stError_ppm <- NA

herbariumData31_78 <- pivot_longer(herbariumData31_78, cols=c(3:17), names_to="metabolite", values_to="concentration_ppm")
herbariumData31_78$stError_ppm <- NA



write.csv(herbariumData1_30, file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200214_herbarium1_30.csv")
write.csv(herbariumData31_78, file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200812_herbarium31_78.csv")