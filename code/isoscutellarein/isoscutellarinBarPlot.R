library(tidyverse)
library(ggrepel)
library(cowplot)

# Load data from .csv files
isoscutellarin <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20190813_isoscutellarin.csv")

# Combine all data into a single data frame
allData <- isoscutellarin

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa_071119", "racemosa_MS", "racemosa_SC", "hastifolia", "havanensis", "arenicola", "suffrutescens"), collapse="|")
excludeOrgans <- paste(c("flowers"), collapse="|")
excludeMetabolites <- paste(c("acteoside"), collapse="|")
allData <- allData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) %>%
  filter(!grepl(excludeMetabolites, metabolite))

# Fix naming errors
allData$species[allData$species=="racemosa_RNAseq"] <- "racemosa"
allData$species[allData$species=="leonardii"] <- "parvula"

# Change classifiers (species, organs, metabolites) into factors
allData$species <- factor(allData$species)
allData$organ <- factor(allData$organ)
allData$metabolite <- factor(allData$metabolite)

allData$metabolite <- as.character(allData$metabolite)
allData$metabolite[allData$metabolite=="isoscutellarin"] <- "isoscutellarein 8-G"
allData$metabolite <- factor(allData$metabolite)

# Capitalize first letter of each flavonoid name
capString <- function(string) {
  c <- strsplit(string, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
}
allData$metabolite <- as.character(allData$metabolite)
allData$metabolite <- sapply(allData$metabolite, capString)

# Set order of organs to appear in plot
allData$organ <- factor(allData$organ, levels=c("leaves", "stems", "roots"))


allData$signif <- NA
allData$signif_Y <- NA
for(row in 1:nrow(allData)){
  if(allData$species[row] != "baicalensis"){
    baiSE <- allData$stError_peakArea[allData$species == "baicalensis" & allData$metabolite == allData$metabolite[row] & allData$organ == allData$organ[row]]
    baiConc <- allData$peakArea[allData$species == "baicalensis" & allData$metabolite == allData$metabolite[row] & allData$organ == allData$organ[row]]
    SEsum <- sqrt(allData$stError_peakArea[row]^2 + baiSE^2) #https://www.graphpad.com/support/faq/the-standard-error-of-the-difference-between-two-means/
    z <- as.numeric(abs(allData$peakArea[row] - baiConc)/SEsum)
    pval <- exp(-0.717*z - 0.416*z^2) #https://www.bmj.com/content/343/bmj.d2304
    if(pval < 0.05 & !is.nan(pval)){
      allData$signif[row] <- "*"
      nudge_Y <- max(allData$peakArea[allData$metabolite == allData$metabolite[row]])*0.025
      allData$signif_Y[row] <- allData$peakArea[row] + allData$stError_peakArea[row] + nudge_Y
    }
  }
}

allData$stError_peakArea[allData$peakArea == 0] <- NA

createIndividualBars <- function(allData, indMetabolite, axisLabels=TRUE, legend=TRUE){
  graphData <- allData %>%
    filter(metabolite==indMetabolite)
  
  graphData$species <- paste("S.", graphData$species)
  graphData$species <- str_wrap(graphData$species, width=16)
  graphData$species <- factor(graphData$species, levels=c("S. baicalensis", "S. altissima", "S. barbata", "S. parvula", "S. racemosa", "S. tournefortii", "S. wrightii"))
  
  indBarPlot <- ggplot(data=graphData, mapping=aes(x=species, y=peakArea, fill=organ)) +
    geom_col(position="dodge") +
    geom_errorbar(mapping=aes(ymin=peakArea-stError_peakArea, ymax=peakArea+stError_peakArea), color="black", width=0.35, size=0.35, position=position_dodge(0.9)) +
    geom_text(mapping=aes(y=signif_Y, label=signif), position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#47acff", "#62c44d", "#ff8d4f"), name="Organ:", breaks=c("leaves", "stems", "roots"), labels=c("Leaf     ", "Stem     ", "Root")) +
    labs(y=paste(indMetabolite, "[peak area]")) +
    theme_classic() +
    theme(panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
          axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text.y=element_text(size=8, color="black"),
          axis.text.x=element_text(size=8, face="italic", color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(4, 0, 0, 0)),
          legend.position="top", legend.title=element_text(size=8, face="bold"), legend.text=element_text(size=8), legend.key.size=unit(0.35, "cm"), legend.box.margin=margin(-5,0,-8,0),
          axis.ticks=element_line(color="black")
    )
  
  if(axisLabels == FALSE){
    indBarPlot <- indBarPlot + theme(axis.text.x=element_blank())
  }
  
  if(legend == FALSE){
    indBarPlot <- indBarPlot + theme(legend.position="none")
  }
  
  
  return(indBarPlot)
}

isoscutellarinPlot <- createIndividualBars(allData, "Isoscutellarein 8-G", axisLabels=TRUE, legend=TRUE)
ggsave(filename="C:/Users/Bryce/Research/scutellariaMetabolites/figures/0-isoscutellarein/Figure_6.pdf",
       plot=isoscutellarinPlot,
       device=pdf(),
       width=4, height=3, units="in")