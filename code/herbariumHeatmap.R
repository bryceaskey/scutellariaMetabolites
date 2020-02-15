# Method to create heatmap from herbarium sample HPLC data
library(ggplot2)
library(plyr)
library(dplyr)
library(tibble)
library(cowplot)
library(ggforce)
library(stats)
library(tidyr)

# Read metabolite data from .csv file -------------------------------------------------------------
rawData <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/metaboliteDataHerbarium1-30.csv", na.strings=c("", "NA"), header=TRUE)
rawData[, 1] <- as.character(rawData[, 1])

# Define functions for interpreting injection names -----------------------------------------------
herbariumList <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/herbariumList.csv")

getSampleName <- function(injectionName, nameList){
  sampleName <- nameList[injectionName==nameList$Label, 3]
  return(as.character(sampleName))
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

rawData <- rawData[!grepl("15mix", rawData$injection), ]
row.names(rawData) <- 1:nrow(rawData)

# Create data frame with processed data -----------------------------------------------------------
allData <- data.frame(
  variety=sapply(rawData$injection[29:72], getSampleName, nameList=herbariumList),
  apigenin=sapply(rawData$apigenin[29:72], ppmConversion, metaboliteName="apigenin", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  apigeninG=sapply(rawData$apigeninG[29:72], ppmConversion, metaboliteName="apigeninG", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  scutellarein=sapply(rawData$scutellarein[29:72], ppmConversion, metaboliteName="scutellarein", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  scutellarin=sapply(rawData$scutellarin[29:72], ppmConversion, metaboliteName="scutellarin", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  hispidulin=sapply(rawData$hispidulin[29:72], ppmConversion, metaboliteName="hispidulin", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  hispidulinG=sapply(rawData$hispidulinG[29:72], ppmConversion, metaboliteName="hispidulinG", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  chrysin=sapply(rawData$chrysin[29:72], ppmConversion, metaboliteName="chrysin", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  chrysinG=sapply(rawData$chrysinG[29:72], ppmConversion, metaboliteName="chrysinG", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  wogonin=sapply(rawData$wogonin[29:72], ppmConversion, metaboliteName="wogonin", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  wogonoside=sapply(rawData$wogonoside[29:72], ppmConversion, metaboliteName="wogonoside", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  baicalein=sapply(rawData$baicalein[29:72], ppmConversion, metaboliteName="baicalein", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  baicalin=sapply(rawData$baicalin[29:72], ppmConversion, metaboliteName="baicalin", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  oroxylinA=sapply(rawData$oroxylinA[29:72], ppmConversion, metaboliteName="oroxylinA", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  oroxyloside=sapply(rawData$oroxyloside[29:72], ppmConversion, metaboliteName="oroxyloside", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData),
  acetoside=sapply(rawData$acetoside[29:72], ppmConversion, metaboliteName="acetoside", mix1Rows=11:18, mix2Rows=20:27, rawData=rawData)
)

# Calculate means for herbarium varieties with >1 samples  ----------------------------------------
allData <- aggregate(allData[ , c(2:16)], by=list(allData$variety), mean, na.rm=TRUE)
colnames(allData)[1] <- "variety"

# Convert data into tidy format for for plotting with ggplot --------------------------------------
allData <- gather(allData, "metabolite", "concentration", 2:16)

# Reorder metabolites as they will appear in heatmap ----------------------------------------------
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein", "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")
allData$metabolite <- factor(allData$metabolite, levels=metaboliteOrder)

# Create heatmap ----------------------------------------------------------------------------------
heatmap <- ggplot(data=allData) +
  geom_raster(mapping=aes(x=variety, y=metabolite, fill=concentration)) +
  scale_fill_gradientn(colours=c("#FFFFFFFF", "#0066CC")) +
  coord_fixed() +
  labs(title="Herbarium sample 1-30 flavonoid concetrations", x="Variety", y="Metabolite", fill="Concentration (ppm)") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), text=element_text(size=16))
return(heatmap)
print(heatmap)

