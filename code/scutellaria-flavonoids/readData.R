# Method to prepare and format metabolite data from raw .csv file
library(plyr)
library(dplyr)
library(tibble)

# Read metabolite data from .csv file ----
rawData <- read.csv(file="...", header=TRUE)
rawData[, 1] <- as.character(rawData[, 1])

# Define functions for interpreting injection names ----
abbrevNames <- data.frame(
  abbrev=c("HV", "AC", "AS", "BT", "TT", "HF", "BL", "LD", "RMSEQ", "R071119", "R_MS", "R_SC"),
  fullName=c("havenesis", "arenicola", "altissima", "barbata", "tournefortii", "hastafolia", "baicalensis", "leonardii", "RNA Seq", "racemosa 071119", "racemosa MS", "racemosa SC")
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

# Define function to convert peak area data into ppm using calibration samples ----
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

# Create data frame with processed data ----
allData <- data.frame(
  variety=sapply(rawData$injectionName[23:112], getSampleName),
  replicate=sapply(rawData$injectionName[23:112], getSampleRep),
  organ=sapply(rawData$injectionName[23:112], getSampleOrgan),
  apigenin=sapply(rawData$apigenin[23:112], ppmConversion, metaboliteName="apigenin", rawData=rawData),
  apigeninG=sapply(rawData$apigeninG[23:112], ppmConversion, metaboliteName="apigeninG", rawData=rawData),
  scutellarein=sapply(rawData$scutellarein[23:112], ppmConversion, metaboliteName="scutellarein", rawData=rawData),
  scutellarin=sapply(rawData$scutellarin[23:112], ppmConversion, metaboliteName="scutellarin", rawData=rawData),
  hispidulin=sapply(rawData$hispidulin[23:112], ppmConversion, metaboliteName="hispidulin", rawData=rawData),
  hispidulinG=sapply(rawData$hispidulinG[23:112], ppmConversion, metaboliteName="hispidulinG", rawData=rawData),
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

# Calculate mean and standard error ----
# Make one data frame for each metabolite, and combine into list.
# Columns: variety / organ / meanConc / stError / upperSD / lowerSD
listData <- list()
for(metName in colnames(allData)[4:ncol(allData)]){
  df <- ddply(allData, c("variety", "organ"), summarise, meanConc=mean(get(metName)), stError=sd(get(metName))/sqrt(length(get(metName))))
  df <- transform(df, lowerSD=meanConc-stError, upperSD=meanConc+stError)
  assign(metName, df)
}
listData <- list(apigenin, apigeninG, scutellarein, scutellarin, hispidulin, hispidulinG, chrysin, chrysinG, wogonin, wogonoside, baicalein, baicalin, oroxylinA, oroxyloside, acetoside)
names(listData) <- colnames(allData)[4:ncol(allData)]
rm(df, list=names(listData))

# Data structure for creating figures with in ggplot ----
allData <- do.call(rbind, listData)
allData$metabolite <- factor(names(listData)[rep(1:length(listData), each=sapply(listData, nrow)[1])])
rownames(allData) <- seq(1, nrow(allData))

# Remove racemosa and barbata data ----
allData <- allData %>%
  filter(!grepl("racemosa MS", variety)) %>%
  filter(!grepl("racemosa SC", variety)) %>%
  filter(!grepl("racemosa 071119", variety)) %>%
  filter(!grepl("RNA Seq", variety))
allData$variety <- factor(allData$variety)

# Adjust allData structure for easier plotting ----
varietyOrder <- c("altissima", "arenicola", "baicalensis", "barbata", "hastafolia", "havenesis", "leonardii", "tournefortii", "racemosa 071119", "racemosa MS", "racemosa SC", "RNA Seq")
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein", "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")
organOrder <- c("Roots", "Shoots", "Leaves")
allData$variety <- factor(allData$variety, levels=varietyOrder)
colnames(allData)[1] <- "species"
allData$metabolite <- factor(allData$metabolite, levels=metaboliteOrder)
allData$organ <- factor(allData$organ, levels=organOrder)
allData$metNum <- as.numeric(allData$metabolite) #create new column w/ factor nums - for bar plot section labeling
allData$species <- factor(allData$species)

# Fix error in naming ----
levels(allData$species)[levels(allData$species)=="hastafolia"] <- "hastifolia"

# Order according to metabolite number ----
allData <- allData[order(allData$metNum), ]