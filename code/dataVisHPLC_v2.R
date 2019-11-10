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
rawData <- read.csv(file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/metaboliteData.csv", header=TRUE)
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

# TODO: Reconsider data structure. Make one data frame for each metabolite, and combine into list?
# Columns: meanConc / stError / upperSD / lowerSD

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

# Better data structure to make scaled point figure with ------------------------------------------
allData <- do.call(rbind, listData)
allData$metabolite <- factor(names(listData)[rep(1:length(listData), each=sapply(listData, nrow)[1])])
rownames(allData) <- seq(1, nrow(allData))

# Data structure for creating organ-specific raster plots -----------------------------------------
organData <- list(
  Leaves=subset(allData, organ=="Leaves", select=c("variety", "meanConc", "metabolite")),
  Shoots=subset(allData, organ=="Shoots", select=c("variety", "meanConc", "metabolite")),
  Roots=subset(allData, organ=="Roots", select=c("variety", "meanConc", "metabolite"))
)

# Scale mean values for each metabolite as a % of max value ---------------------------------------
# Define function to work on dataframe for a single organ
normalizeValues <- function(df){
  normData <- data.frame()
  for(met in levels(df$metabolite)){
    metData <- subset(df, metabolite==met)
    maxConc <- max(metData$meanConc)
    metData <- transform(metData, normMean=meanConc/maxConc)
    normData <- rbind(normData, metData)
    
  }
  return(normData)
}

# Apply function to each organ
organData <- list(
  Leaves=normalizeValues(organData$Leaves),
  Shoots=normalizeValues(organData$Shoots),
  Roots=normalizeValues(organData$Roots)
)

# Create raster plot (i.e. heatmap) for each organ ------------------------------------------------------------------------
# Remove Havenesis, move Tourmetii next to Lateriflora, and adjust metabolite order to match pathway
varietyOrder <- c("Altissima", "Arenicola", "Baicalensis", "Barbata", "Hastafolia", "Lateriflora", "Tourmetii", "Racemosa 071119", "Racemosa MS", "Racemosa SC", "RNA Seq")
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein", "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")
organData$Leaves <- subset(organData$Leaves, variety!="Havenesis")
organData$Leaves$variety <- factor(organData$Leaves$variety, levels=varietyOrder)
organData$Leaves$metabolite <- factor(organData$Leaves$metabolite, levels=metaboliteOrder)
organData$Shoots <- subset(organData$Shoots, variety!="Havenesis")
organData$Shoots$variety <- factor(organData$Shoots$variety, levels=varietyOrder)
organData$Shoots$metabolite <- factor(organData$Shoots$metabolite, levels=metaboliteOrder)
organData$Roots <- subset(organData$Roots, variety!="Havenesis")
organData$Roots$variety <- factor(organData$Roots$variety, levels=varietyOrder)
organData$Roots$metabolite <- factor(organData$Roots$metabolite, levels=metaboliteOrder)

rootHeatmap <- ggplot(data=organData$Roots, mapping=aes(x=variety, y=metabolite, fill=meanConc)) +
  geom_raster() +
  scale_fill_gradientn(colours=c("#FFFFFFFF","#FF3333")) +
  coord_fixed() +
  labs(title="Root metabolites", x="Variety", y="Metabolite", fill="Conc (ppm)") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

shootHeatmap <- ggplot(data=organData$Shoots, mapping=aes(x=variety, y=metabolite, fill=meanConc)) +
  geom_raster() +
  scale_fill_gradientn(colours=c("#FFFFFFFF","#009900")) +
  coord_fixed() + 
  labs(title="Shoot metabolites", x="Variety", y="Metabolite", fill="Conc (ppm)") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

leafHeatmap <- ggplot(data=organData$Leaves, mapping=aes(x=variety, y=metabolite, fill=meanConc)) +
  geom_raster() +
  scale_fill_gradientn(colours=c("#FFFFFFFF","#0066CC")) +
  coord_fixed() +
  labs(title="Leaf metabolites", x="Variety", y="Metabolite", fill="Conc (ppm)") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

plot_grid(rootHeatmap, shootHeatmap, leafHeatmap, nrow=1, ncol=3)

# Method to create scaled pie charts --------------------------------------------------------------
altissima <- ggplot(transform(subset(allData, variety=="Altissima"), organ=factor(organ, levels=c("Roots", "Shoots", "Leaves")))) +
  geom_bar(mapping=aes(x="", y=meanConc, fill=metabolite), width=1, stat="identity", position="fill", color="white") +
  coord_polar("y", start=0) +
  facet_grid(~organ) +
  theme(legend.position="bottom", legend.direction="horizontal")
# Remove rows with meanConc = 0, but maintain consistent color scheme across all plots
# Colors must be first be manually assigned to metabolites
metaboliteColors <- c("#8B0000", "#DC143C", "#FF7F50", "#FFD700", "#B8860B", "#BDB76B", "#808000", "#9ACD32", "#2E8B57", "#66CDAA", "#2F4F4F", "#008080", "#4682B4", "#8A2BE2", "#8B008B")
names(metaboliteColors) <- metaboliteOrder

altissimaRootData <- subset(subset(subset(allData, variety=="Altissima"), organ=="Roots"), meanConc>0)
altissimaShootData <- subset(subset(subset(allData, variety=="Altissima"), organ=="Shoots"), meanConc>0)
altissimaLeafData <- subset(subset(subset(allData, variety=="Altissima"), organ=="Leaves"), meanConc>0)

# Calculate position of labels in pie charts ------------------------------------------------------
altissimaRootData <- altissimaRootData %>%
  mutate(end=2*pi*cumsum(meanConc)/sum(meanConc),
   start=lag(end, default=0),
   middle=0.5*(start+end),
   hjust=ifelse(middle>pi, 1, 0),
   vjust=ifelse(middle<pi/2 | middle>3*pi/2, 0, 1))

altissimaRoot <- ggplot(altissimaRootData) +
  geom_arc_bar(mapping=aes(x0=0, y0=0, r0=0, r=1, start=start, end=end, fill=metabolite), show.legend=FALSE, color="white") +
  geom_text(mapping=aes(x=1.05*sin(middle), y=1.05*cos(middle), label=metabolite, hjust=hjust, vjust=vjust)) +
  coord_fixed() +
  scale_x_continuous(limits=c(-1.5, 1.5), name="", breaks=NULL, labels=NULL) +
  scale_y_continuous(limits=c(-1.5, 1.5), name="", breaks=NULL, labels=NULL) +
  scale_fill_manual(values=metaboliteColors) +
  theme(panel.background=element_blank(), plot.background=element_blank(), axis.line=element_blank())

altissimaShootData <- altissimaShootData %>%
  mutate(end=2*pi*cumsum(meanConc)/sum(meanConc),
   start=lag(end, default=0),
   middle=0.5*(start+end),
   hjust=ifelse(middle>pi, 1, 0),
   vjust=ifelse(middle<pi/2 | middle>3*pi/2, 0, 1))

altissimaShoot <- ggplot(altissimaShootData) +
  geom_arc_bar(mapping=aes(x0=0, y0=0, r0=0, r=1, start=start, end=end, fill=metabolite), show.legend=FALSE, color="white") +
  geom_text(mapping=aes(x=1.05*sin(middle), y=1.05*cos(middle), label=metabolite, hjust=hjust, vjust=vjust)) +
  coord_fixed() +
  scale_x_continuous(limits=c(-1.5, 1.5), name="", breaks=NULL, labels=NULL) +
  scale_y_continuous(limits=c(-1.5, 1.5), name="", breaks=NULL, labels=NULL) +
  scale_fill_manual(values=metaboliteColors) +
  theme(panel.background=element_blank(), plot.background=element_blank(), axis.line=element_blank())

altissimaLeafData <- altissimaLeafData %>%
  mutate(end=2*pi*cumsum(meanConc)/sum(meanConc),
   start=lag(end, default=0),
   middle=0.5*(start+end),
   hjust=ifelse(middle>pi, 1, 0),
   vjust=ifelse(middle<pi/2 | middle>3*pi/2, 0, 1))

altissimaLeaf <- ggplot(altissimaLeafData) +
  geom_arc_bar(mapping=aes(x0=0, y0=0, r0=0, r=1, start=start, end=end, fill=metabolite), show.legend=FALSE, color="white") +
  geom_text(mapping=aes(x=1.05*sin(middle), y=1.05*cos(middle), label=metabolite, hjust=hjust, vjust=vjust)) +
  coord_fixed() +
  scale_x_continuous(limits=c(-1.5, 1.5), name="", breaks=NULL, labels=NULL) +
  scale_y_continuous(limits=c(-1.5, 1.5), name="", breaks=NULL, labels=NULL) +
  scale_fill_manual(values=metaboliteColors) +
  theme(panel.background=element_blank(), plot.background=element_blank(), axis.line=element_blank())

# TODO: Make position calculation and pie chart generation into functions to simplify large-scale implementation.
# TODO: Replace metabolite names with numbers in pie chart, and add corresponding number to label in legend.

legend <- cowplot::get_legend(altissima)

pieSizes <- ddply(allData, c("variety", "organ"), summarise, area=sum(meanConc))
pieSizes <- transform(pieSizes, radius=sqrt(area/pi))

print(plot_grid(altissimaRoot, altissimaShoot, altissimaLeaf, 
  nrow=1, ncol=3,
  labels=c("Root", "Shoot", "Leaf"), label_x=0, label_y=1, label_size=14,
  rel_widths=c(pieSizes$radius[2], pieSizes$radius[3], pieSizes$radius[1]),
  rel_heights=c(pieSizes$radius[2], pieSizes$radius[3], pieSizes$radius[1])
))
