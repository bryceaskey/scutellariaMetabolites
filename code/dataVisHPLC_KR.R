# Injection names should be standardized to match the following format:
# Variety Abbreviation - Replicate Number - Plant Organ
# e.g. BL-1-R = baicalensis - 1st replicate - Root

library(ggplot2)
library(plyr)
library(dplyr)
library(tibble)
library(cowplot)
library(ggforce)

# Read metabolite data from .csv file -------------------------------------------------------------
rawData <- read.csv(file="C:/Users/bca08_000/Documents/scutellariaMetabolites/data/metaboliteData.csv", header=TRUE)
rawData[, 1] <- as.character(rawData[, 1])

rawDataKR <- read.csv(file="C:/Users/bca08_000/Documents/scutellariaMetabolites/data/2020-01-17_metaboliteDataKR.csv", header=TRUE)
rawDataKR [, 1] <- as.character(rawDataKR[, 1])
colnames(rawDataKR)[1] <- "injectionName"

# Define functions for interpreting injection names -----------------------------------------------
abbrevNames <- data.frame(
  abbrev=c("HV", "AC", "AS", "BT", "TT", "HF", "BL", "LD", "RMSEQ", "R071119", "R_MS", "R_SC"),
  fullName=c("Havenesis", "Arenicola", "Altissima", "Barbata", "Tournefortii", "Hastifolia", "Baicalensis", "Leonardii", "RNA Seq", "Racemosa 071119", "Racemosa MS", "Racemosa SC")
)

abbrevNamesKR <- data.frame(
  abbrev=c("IC", "BT", "DE", "IN", "PA", "ST"),
  fullName=c("KR_Indica", "KR_Barbata", "KR_Dependens", "KR_Insignis", "KR_Pekinesis", "KR_Strigillosa")
)

getSampleName <- function(injectionName, abbrevs){
  sampleAbbrev <- strsplit(injectionName, "-")[[1]][1]
  sampleName <- abbrevs[grep(sampleAbbrev, abbrevs$abbrev), 2]
  return(sampleName)
}

getSampleRep <- function(injectionName){
  sampleRep <- strsplit(injectionName, "-")[[1]][2]
  return(sampleRep)
}

abbrevOrgans <- data.frame(
  abbrev=c("L", "S", "R", "F"),
  fullOrgan=c("Leaves", "Shoots", "Roots", "Flowers")
)

getSampleOrgan <- function(injectionName, abbrevs){
  organAbbrev <- strsplit(injectionName, "-")[[1]][3]
  sampleOrgan <- abbrevOrgans[grep(organAbbrev, abbrevOrgans$abbrev), 2]
  return(sampleOrgan)
}

# Define function to convert peak area data into ppm using calibration samples --------------------
ppmConversion <- function(peakArea, metaboliteName, mix1Rows, mix2Rows, rawData=rawData){
  if(sum(grepl(metaboliteName, colnames(rawData)[2:9]))==1){ #metabolite is in mix 1
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][mix1Rows])
  }else{ #metabolite is in mix 2
    calibrations <- data.frame(ppm=c(0.1, 0.5, 1, 5, 10, 25, 50, 100), area=rawData[[metaboliteName]][mix2Rows])
  }
  convFactor <- lm(formula=ppm~area+0, data=calibrations)[[1]][[1]]
  ppm <- peakArea*convFactor
  return(ppm)
}

# Create data frame with processed data -----------------------------------------------------------
allData <- data.frame(
  variety=sapply(rawData$injectionName[23:112], getSampleName, abbrevs=abbrevNames),
  replicate=sapply(rawData$injectionName[23:112], getSampleRep),
  organ=sapply(rawData$injectionName[23:112], getSampleOrgan, abbrevs=abbrevOrgans),
  apigenin=sapply(rawData$apigenin[23:112], ppmConversion, metaboliteName="apigenin", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  apigeninG=sapply(rawData$apigeninG[23:112], ppmConversion, metaboliteName="apigeninG", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  scutellarein=sapply(rawData$scutellarein[23:112], ppmConversion, metaboliteName="scutellarein", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  scutellarin=sapply(rawData$scutellarin[23:112], ppmConversion, metaboliteName="scutellarin", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  hispidulin=sapply(rawData$hispidulin[23:112], ppmConversion, metaboliteName="hispidulin", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  hispidulinG=sapply(rawData$hispidulinG[23:112], ppmConversion, metaboliteName="hispidulinG", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  chrysin=sapply(rawData$chrysin[23:112], ppmConversion, metaboliteName="chrysin", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  chrysinG=sapply(rawData$chrysinG[23:112], ppmConversion, metaboliteName="chrysinG", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  wogonin=sapply(rawData$wogonin[23:112], ppmConversion, metaboliteName="wogonin", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  wogonoside=sapply(rawData$wogonoside[23:112], ppmConversion, metaboliteName="wogonoside", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  baicalein=sapply(rawData$baicalein[23:112], ppmConversion, metaboliteName="baicalein", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  baicalin=sapply(rawData$baicalin[23:112], ppmConversion, metaboliteName="baicalin", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  oroxylinA=sapply(rawData$oroxylinA[23:112], ppmConversion, metaboliteName="oroxylinA", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  oroxyloside=sapply(rawData$oroxyloside[23:112], ppmConversion, metaboliteName="oroxyloside", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData),
  acetoside=sapply(rawData$acetoside[23:112], ppmConversion, metaboliteName="acetoside", mix1Rows=5:12, mix2Rows=14:21, rawData=rawData)
)

rawDataKR <- rawDataKR[!grepl("15mix", rawDataKR$injectionName), ]
row.names(rawDataKR) <- 1:nrow(rawDataKR)

allDataKR <- data.frame(
  variety=sapply(rawDataKR$injectionName[29:88], getSampleName, abbrevs=abbrevNamesKR),
  replicate=sapply(rawDataKR$injectionName[29:88], getSampleRep),
  organ=sapply(rawDataKR$injectionName[29:88], getSampleOrgan, abbrevs=abbrevOrgans),
  apigenin=rawDataKR$apigenin[29:88],
  apigeninG=rawDataKR$apigeninG[29:88],
  scutellarein=rawDataKR$scutellarein[29:88],
  scutellarin=rawDataKR$scutellarin[29:88],
  hispidulin=rawDataKR$hispidulin[29:88],
  hispidulinG=rawDataKR$hispidulinG[29:88],
  chrysin=rawDataKR$chrysin[29:88],
  chrysinG=rawDataKR$chrysinG[29:88],
  wogonin=rawDataKR$wogonin[29:88],
  wogonoside=rawDataKR$wogonoside[29:88],
  baicalein=rawDataKR$baicalein[29:88],
  baicalin=rawDataKR$baicalin[29:88],
  oroxylinA=rawDataKR$oroxylinA[29:88],
  oroxyloside=rawDataKR$oroxyloside[29:88],
  acetoside=rawDataKR$acetoside[29:88]
)

# Calculate mean and standard error -------------------------------------------------------------
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

listDataKR <- list()
for(metName in colnames(allDataKR)[4:ncol(allDataKR)]){
  df <- ddply(allDataKR, c("variety", "organ"), summarise, meanConc=mean(get(metName)), stError=sd(get(metName))/sqrt(length(get(metName))))
  df <- transform(df, lowerSD=meanConc-stError, upperSD=meanConc+stError)
  assign(metName, df)
}
listDataKR <- list(apigenin, apigeninG, scutellarein, scutellarin, hispidulin, hispidulinG, chrysin, chrysinG, wogonin, wogonoside, baicalein, baicalin, oroxylinA, oroxyloside, acetoside)
names(listDataKR) <- colnames(allDataKR)[4:ncol(allDataKR)]
rm(df, list=names(listData))

# Data structure for creating figures with in ggplot ----------------------------------------------
allData <- do.call(rbind, listData)
allData$metabolite <- factor(names(listData)[rep(1:length(listData), each=sapply(listData, nrow)[1])])
rownames(allData) <- seq(1, nrow(allData))

allDataKR <- do.call(rbind, listDataKR)
allDataKR$metabolite <- factor(names(listDataKR)[rep(1:length(listDataKR), each=sapply(listDataKR, nrow)[1])])
rownames(allDataKR) <- seq(1, nrow(allDataKR))

combinedData <- rbind(allData, allDataKR)

# Adjust combinedData structure for easier plotting ----------------------------------------------------
varietyOrder <- c("Altissima", "Leonardii", "Hastifolia", "Havenesis", "Arenicola", "Tournefortii", "RNA Seq", "Baicalensis", "Barbata", "Racemosa 071119", "Racemosa MS", "Racemosa SC", "KR_Indica", "KR_Barbata", "KR_Dependens", "KR_Insignis", "KR_Pekinesis", "KR_Strigillosa")
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein", "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")
organOrder <- c("Roots", "Shoots", "Leaves", "Flowers")
combinedData$variety <- factor(combinedData$variety, levels=varietyOrder)
combinedData$metabolite <- factor(combinedData$metabolite, levels=metaboliteOrder)
combinedData$organ <- factor(combinedData$organ, levels=organOrder)
combinedData$metNum <- as.numeric(combinedData$metabolite) #create new column w/ factor nums - for pie chart sector labeling

# Remove non-triplicate Racemosa samples, and keep RNA Seq data as true Racemosa data -------------
combinedData <- filter(combinedData, variety != "Racemosa 071119" & variety != "Racemosa MS" & variety != "Racemosa SC")
combinedData$variety <- factor(combinedData$variety, levels=c(levels(combinedData$variety), "Racemosa"))
combinedData$variety[combinedData$variety=="RNA Seq"] <- "Racemosa"
combinedData$variety <- factor(combinedData$variety, levels=c("Altissima", "Leonardii", "Hastifolia", "Havenesis", "Arenicola", "Tournefortii", "Baicalensis", "Barbata", "Racemosa", "KR_Indica", "KR_Barbata", "KR_Dependens", "KR_Insignis", "KR_Pekinesis", "KR_Strigillosa"))

# Define function to create raster plot (i.e. heatmap) for each organ -----------------------------
createHeatmap <- function(combinedData, plantOrgan){
  if(plantOrgan=="Roots"){
    colorScale <- c("#FFFFFFFF", "#FF3333")
    plotTitle <- "Root metabolites"
  }else if(plantOrgan=="Shoots"){
    colorScale <- c("#FFFFFFFF", "#009900")
    plotTitle <- "Shoot metabolites"
  }else if(plantOrgan=="Leaves"){
    colorScale <- c("#FFFFFFFF", "#0066CC")
    plotTitle <- "Leaf metabolites"
  }else{
    colorScale <- c("#FFFFFFFF", "#FCBA03")
    plotTitle <- "Flower metabolites"
  }
  organData <- filter(combinedData, organ==plantOrgan)
  heatmap <- ggplot(data=organData) +
    geom_raster(mapping=aes(x=variety, y=metabolite, fill=meanConc)) +
    scale_fill_gradientn(colours=colorScale) +
    coord_fixed() +
    labs(title=plotTitle, x="Variety", y="Metabolite", fill="Conc (ppm)") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  return(heatmap)
}

print(plot_grid(
  createHeatmap(combinedData, "Leaves"),
  createHeatmap(combinedData, "Shoots"),
  createHeatmap(combinedData, "Roots"),
  nrow=1, ncol=3,
  rel_widths=c(1, 1, 1), rel_heights=c(1, 1, 1)))