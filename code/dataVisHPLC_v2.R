# Injection names should be standardized to match the following format:
# Variety Abbreviation - Replicate Number - Plant Organ
# e.g. BL-1-R = baicalensis - 1st replicate - Root

library(ggplot2)
library(plyr)
library(dplyr)
library(tibble)
library(cowplot)
library(ggforce)
library(grid)
library(gridExtra)

# Read metabolite data from .csv file -------------------------------------------------------------
rawData <- read.csv(file="C:/Users/bca08_000/Documents/scutellariaMetabolites/data/metaboliteData.csv", header=TRUE)
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

# Data structure for creating figures with in ggplot ----------------------------------------------
allData <- do.call(rbind, listData)
allData$metabolite <- factor(names(listData)[rep(1:length(listData), each=sapply(listData, nrow)[1])])
rownames(allData) <- seq(1, nrow(allData))

# Adjust allData structure for easier plotting ----------------------------------------------------
varietyOrder <- c("Altissima", "Arenicola", "Baicalensis", "Barbata", "Hastafolia", "Lateriflora", "Tourmetii", "Racemosa 071119", "Racemosa MS", "Racemosa SC", "RNA Seq")
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin", "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein", "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")
organOrder <- c("Roots", "Shoots", "Leaves")
allData <- filter(allData, variety!="Havenesis") #remove Havenesis
allData$variety <- factor(allData$variety, levels=varietyOrder)
allData$metabolite <- factor(allData$metabolite, levels=metaboliteOrder)
allData$organ <- factor(allData$organ, levels=organOrder)
allData$metNum <- as.numeric(allData$metabolite) #create new column w/ factor nums - for pie chart sector labeling

# Define function to create raster plot (i.e. heatmap) for each organ -----------------------------
createHeatmap <- function(allData, plantOrgan){
  if(plantOrgan=="Roots"){
    colorScale <- c("#FFFFFFFF", "#FF3333")
    plotTitle <- "Root metabolites"
  }else if(plantOrgan=="Shoots"){
    colorScale <- c("#FFFFFFFF", "#009900")
    plotTitle <- "Shoot metabolites"
  }else{
    colorScale <- c("#FFFFFFFF", "#0066CC")
    plotTitle <- "Leaf metabolites"
  }
  organData <- filter(allData, organ==plantOrgan)
  heatmap <- ggplot(data=organData) +
    geom_raster(mapping=aes(x=variety, y=metabolite, fill=meanConc)) +
    scale_fill_gradientn(colours=colorScale) +
    coord_fixed() +
    labs(title=plotTitle, x="Variety", y="Metabolite", fill="Conc (ppm)") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  return(heatmap)
}


# Set colors to be used for metabolites across all plots ------------------------------------------
metaboliteColors <- c("#8B0000", "#DC143C", "#FF7F50", "#FFD700", "#B8860B", "#BDB76B", "#808000", "#9ACD32", "#2E8B57", "#66CDAA", "#2F4F4F", "#008080", "#4682B4", "#8A2BE2", "#8B008B")
names(metaboliteColors) <- metaboliteOrder

# Define function to create color legend for scaled pie charts ------------------------------------
createLegend <- function(allData, metaboliteColors, legendOrientation="horizontal"){
  x <- ggplot(filter(allData, variety==levels(allData$variety)[1])) +
    geom_bar(mapping=aes(x="", y=meanConc, fill=metabolite), stat="identity") +
    theme(legend.position="bottom", legend.direction=legendOrientation, legend.text=element_text(size=10)) +
    scale_fill_manual(values=metaboliteColors, labels=paste(seq(1, 15), ". ", names(metaboliteColors), sep=""))
  legend <- get_legend(x)
  return(legend)
}

# Define function to create pie charts ------------------------------------------------------------
createPieChart <- function(allData, metaboliteColors, plantVariety, plantOrgan, size){
  subData <- filter(allData, variety==plantVariety & organ==plantOrgan & meanConc>0)
  subData <- subData[order(subData$metNum), ]
  subData <- subData %>% mutate(end=2*pi*cumsum(meanConc)/sum(meanConc),
    start=lag(end, default=0),
    middle=0.5*(start+end),
    hjust=ifelse(middle>pi, 1, 0),
    vjust=ifelse(middle<pi | middle>3*pi/2, 0, 1))
  pieChart <- ggplot(data=subData) +
    geom_arc_bar(mapping=aes(x0=0, y0=0, r0=0, r=size, start=start, end=end, fill=metabolite), show.legend=FALSE, color="white") +
    #TODO: Add condition to detect slivers of pie that are directly next to each other - test with geom_text_repel
    geom_text(mapping=aes(x=(size+0.3)*sin(middle), y=(size+0.3)*cos(middle), label=metNum), size=3) +
    geom_segment(mapping=aes(x=(size+0.05)*sin(middle), y=(size+0.05)*cos(middle), xend=(size+0.15)*sin(middle), yend=(size+0.15)*cos(middle)), color="black", size=0.5) +
    coord_fixed() +
    scale_x_continuous(limits=c(-1.325, 1.325), name="", breaks=NULL, labels=NULL) +
    scale_y_continuous(limits=c(-1.325, 1.325), name="", breaks=NULL, labels=NULL) +
    scale_fill_manual(values=metaboliteColors) +
    theme(panel.background=element_blank(), plot.background=element_blank(), plot.margin=unit(c(0, 0, 0, 0), "cm"), axis.line=element_blank())
  return(pieChart)
}

# Define function to create labels for complete figure --------------------------------------------

# TODO: Increase distance of labels from sectors, and add tick marks to middle.

# TODO: If using r value to scale, need to develop method to scale labels with pie charts

# Subset allData to select varieties for plotting with scaled pie charts  -------------------------
plottingData <- filter(allData, variety=="Arenicola" | variety=="Barbata" | variety=="Altissima")
plottingData$variety <- factor(plottingData$variety)

# Calculate size of pies based on total amount of metabolites -------------------------------------
pieSizes <- ddply(plottingData, c("variety", "organ"), summarise, area=sum(meanConc))
pieSizes <- transform(pieSizes, radius=sqrt(area/pi))
pieSizes <- transform(pieSizes, scaledRadius=pieSizes$radius/max(radius))

# For loop to create list of pie charts for plotting ----------------------------------------------
allPies <- vector("list", length=length(levels(plottingData$variety))*3)

i = 0
for(variety in levels(plottingData$variety)){
  for(organ in levels(plottingData$organ)){
    i = i + 1
    allPies[[i]] <- local({
      i <- i
      pieRadius <- pieSizes[which(pieSizes$variety == variety & pieSizes$organ == organ), 5]
      pieChart <- createPieChart(plottingData, metaboliteColors, variety, organ, size=pieRadius)
    })
    names(allPies)[i] <- paste(variety, organ)
  }
}

legend <- createLegend(allData, metaboliteColors)
justPies <- plot_grid(plotlist=allPies, 
  ncol=3)

y.grob <- textGrob("Barbata Arenicola Altissima", 
  gp=gpar(fontface="bold", col="black", fontsize=14), rot=90)
x.grob <- textGrob("Roots Shoots Leaves", 
  gp=gpar(fontface="bold", col="black", fontsize=14))

justPies <- grid.arrange(arrangeGrob(justPies, left=y.grob, top=x.grob))

finalFigure <- plot_grid(justPies, legend,
  nrow=2, ncol=1, 
  rel_heights=c(1, 0.15))

print(finalFigure)