# Create bar graphs for visualization of scutellaria metabolite data.
# Separate bar graph for each metabolite, sorted by variety, then by part of plant (root, shoot, or leaf).
packageList <- c("ggplot2")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0){
  install.packages(newPackages)
}
library(ggplot2)

setwd("C:/Users/bca08_000/Documents/scutellaria")

allData <- read.csv("metaboliteData.csv", header=TRUE, stringsAsFactors=FALSE)[, 1:16]
calibrationData <- read.csv("metaboliteCalibrations.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(calibrationData)[1] <- "ppm"

# prepare data - delete rows with wash and calibration samples, and average triplicates together
metaboliteData <- allData[23:112, ]
sampleName <- vector()
plantOrgan <- vector()
avgMetData <- data.frame()

for(row in seq(1, 81, 3)){
  sampleName <- append(sampleName, strsplit(metaboliteData$injectionName[row], "-")[[1]][1])
  plantOrgan <- append(plantOrgan, strsplit(metaboliteData$injectionName[row], "-")[[1]][3])
  metReplicates <- rbind(metaboliteData[row, 2:16], metaboliteData[row+1, 2:16], metaboliteData[row+2, 2:16])
  avgMetData <- rbind(avgMetData, colMeans(metReplicates))
  #metabolites <- colMeans(subset(metaboliteData, grepl(metaboliteData$injectionName, sampleName[length(sampleName)]), select=2:16))
}
colnames(avgMetData) <- colnames(metaboliteData[, 2:16])

for(row in 82:90){
  sampleName <- append(sampleName, strsplit(metaboliteData$injectionName[row], "-")[[1]][1])
  plantOrgan <- append(plantOrgan, strsplit(metaboliteData$injectionName[row], "-")[[1]][2])
  avgMetData <- rbind(avgMetData, metaboliteData[row, 2:16])
}

metaboliteData <- cbind(sampleName, plantOrgan, avgMetData)

metaboliteData$sampleName <- as.character(metaboliteData$sampleName)
metaboliteData$sampleName <- factor(metaboliteData$sampleName, levels=unique(metaboliteData$sampleName))
metaboliteData$plantOrgan <- as.character(metaboliteData$plantOrgan)
metaboliteData$plantOrgan <- factor(metaboliteData$plantOrgan, levels=unique(metaboliteData$plantOrgan))

for(i in 3:17){
  print(
    ggplot(data=metaboliteData, mapping=aes(x=sampleName, y=get(colnames(metaboliteData)[i]), fill=plantOrgan)) + 
      geom_bar(stat="identity", position="dodge") +
      labs(x="Injection Name", y="Peak Area [mUA*min]", fill="Plant Organ", title=paste("08/13/19", colnames(metaboliteData)[[i]])) +
      theme(plot.title=element_text(hjust=0.5, size=20), axis.title=element_text(size=16), axis.text=element_text(size=12), legend.title=element_text(size=16), legend.text=element_text(size=12)))
}
