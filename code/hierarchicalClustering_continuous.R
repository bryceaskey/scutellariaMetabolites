library(tidyverse)
library(ggplot2)
library(ggdendro)
library(cowplot)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]
herbarium1_30 <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200214_herbarium1_30.csv")[, 2:6]

# Combine all data into a single data frame and change classifiers (species, organs, metabolites)
# into factors
allData <- rbind(fresh, frozenKR, herbarium1_30)
allData$species <- as.factor(allData$species)
allData$organ <- as.factor(allData$organ)
allData$metabolite <- as.factor(allData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastifolia"), collapse = '|')
excludeOrgans <- paste(c("flowers"), collapse = '|')
#excludeMetabolites <- paste(c("chrysinG", "oroxyloside", "baicalin", "wogonoside", "acetoside", "apigeninG", "scutellarin", "hispidulinG"), collapse = '|')
allData <- allData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) #%>%
#filter(!grepl(excludeMetabolites, metabolite))

# Fix naming error
allData$species <- as.character(allData$species)
allData$species[allData$species=="RNA Seq"] <- "racemosa"
allData$species[allData$species=="havenesis"] <- "havanensis"
allData$species[allData$species=="hastafolia"] <- "hastifolia"
allData$species <- as.factor(allData$species)

# Separate organ-specific data from non organ-specific data
organSpData <- filter(allData, organ != "all")
nonOrganSpData <- filter(allData, organ == "all")

# Average together organ-specific data and recombine with non organ-specific data
organSpData_avgd <- organSpData %>%
  group_by(species, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm))
nonOrganSpData$organ <- NULL
heatmapData <- rbind(organSpData_avgd, nonOrganSpData)

# Average together any duplicate data points
heatmapData <- heatmapData %>%
  group_by(species, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm))

# Define function to convert units of ppm to micromol/L
ppm2microM <- function(input_ppm, metaboliteName){
  if(!is.na(input_ppm)){
    if(metaboliteName=="acetoside"){ #PubChem CID: 5281800 
      output_microM <- (input_ppm/624.6)*1000
    }else if(metaboliteName=="apigenin"){ #PubChem CID: 5280443
      output_microM <- (input_ppm/270.24)*1000
    }else if(metaboliteName=="apigeninG"){ #PubChem CID: 5280704
      output_microM <- (input_ppm/432.4)*1000
    }else if(metaboliteName=="baicalein"){ #PubChem CID: 5281605
      output_microM <- (input_ppm/270.24)*1000
    }else if(metaboliteName=="baicalin"){ #PubChem CID: 64982
      output_microM <- (input_ppm/446.4)*1000
    }else if(metaboliteName=="chrysin"){ #PubChem CID: 5281607
      output_microM <- (input_ppm/254.24)*1000
    }else if(metaboliteName=="chrysinG"){ #PubChem CID: 90658886
      output_microM <- (input_ppm/416.4)*1000
    }else if(metaboliteName=="hispidulin"){ #PubChem CID: 5281628
      output_microM <- (input_ppm/300.26)*1000
    }else if(metaboliteName=="hispidulinG"){ #PubChem CID: 5318083
      output_microM <- (input_ppm/462.4)*1000
    }else if(metaboliteName=="oroxylinA"){ #PubChem CID: 5320315
      output_microM <- (input_ppm/284.26)*1000
    }else if(metaboliteName=="oroxyloside"){ #PubChem CID: 14655551
      output_microM <- (input_ppm/460.4)*1000
    }else if(metaboliteName=="scutellarein"){ #PubChem CID: 5281697
      output_microM <- (input_ppm/286.24)*1000
    }else if(metaboliteName=="scutellarin"){ #PubChem CID: 185617
      output_microM <- (input_ppm/462.4)*1000
    }else if(metaboliteName=="wogonin"){ #PubChem CID: 5281703
      output_microM <- (input_ppm/284.26)*1000
    }else if(metaboliteName=="wogonoside"){ #PubChem CID: 3084961
      output_microM <- (input_ppm/460.4)*1000
    }else{
      print(paste("Error: metabolite name", metaboliteName,  "not recognized"))
      output_microM <- NA
    }
  }else{
    output_microM <- NA
  }
  return(output_microM)
}

# Convert concentration and stError from ppm to mM for each data point
concentration_microM <- vector(mode="numeric", length=nrow(heatmapData))
stError_microM <- vector(mode="numeric", length=nrow(heatmapData))
for(i in 1:nrow(heatmapData)){
  concentration_microM[i] <- ppm2microM(heatmapData$concentration_ppm[i], heatmapData$metabolite[i])
  stError_microM[i] <- ppm2microM(heatmapData$stError_ppm[i], heatmapData$metabolite[i])
}
heatmapData$concentration_microM <- concentration_microM
heatmapData$stError_microM <- stError_microM

# Assign metabolite order as determined by heirarchical clustering
metaboliteOrder <- c("chrysin", "chrysinG", "oroxylinA", "oroxyloside", "baicalein",
                     "baicalin", "wogonin", "wogonoside", "acetoside", "apigenin",
                     "apigeninG", "scutellarein", "scutellarin", "hispidulin", "hispidulinG")
heatmapData$metabolite <- factor(heatmapData$metabolite, levels=metaboliteOrder)

# Transform data into wide format to use for heirarchical clustering 
contData <- subset(heatmapData, select=-c(concentration_ppm, stError_ppm, stError_microM))
contData <- contData %>%
  pivot_wider(names_from=metabolite, values_from=concentration_microM) %>%
  remove_rownames %>%
  column_to_rownames(var="species")
contData <- scale(contData)

# Calculate distance matrix
contDist <- dist(contData, method="euclidian")

contCluster <- hclust(d=contDist, method="average")
print(paste("Correlation between cophenetic distances and actual distances:", cor(contDist, cophenetic(contCluster))))
contDenData <- dendro_data(as.dendrogram(contCluster), type="rectangle")

# Load clade data from phylogenetic tree generated from chloroplast genome
cladeData <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/phylo-tree-clades.csv")
speciesList <- vector(mode="character", length=nrow(label(contDenData)))
cladeList <- vector(mode="numeric", length=nrow(label(contDenData)))
for (i in 1:nrow(label(contDenData))){
  denSpecies <- label(contDenData)$label[i]
  speciesList[i] <- denSpecies
  clade <- cladeData$clade[grep(denSpecies, cladeData$species)]
  if (length(clade)>0){
    cladeList[i] <- clade
  }else{
    cladeList[i] <- NA
  }
} 

denCladeData <- data.frame(x=1, y=1:nrow(label(contDenData)), speciesList, cladeList)
denCladeData$cladeList <- factor(denCladeData$cladeList, levels=c(1, 2, 3, 4, 5, 6))

# Create dendrogram
contDenPlot <- ggplot() +
  geom_segment(data=segment(contDenData), mapping=aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(contDenData), mapping=aes(x=x, y=y, label=label, hjust=1), nudge_y=7) +
  coord_flip() +
  scale_y_reverse(expand=c(0.3, 0)) +
  theme_dendro() +
  theme(plot.margin=unit(c(0.00075,0,0.067,0.1),"npc"))

# Adjust species order in heatmap data to match order in dendrogram
speciesOrder <- label(contDenData)$label
heatmapData$species <- factor(heatmapData$species, levels=speciesOrder)

# Create heatmap without species labels
heatmap <- ggplot(data=heatmapData) +
  geom_raster(mapping=aes(x=species, y=metabolite, fill=concentration_microM)) +
  scale_fill_viridis() +
  labs(y="Flavonoid", fill=expression(paste("Conc (", mu, "M)", sep=""))) +
  coord_flip() +
  theme(axis.text.x=element_text(angle=45, hjust=1, color="#000000"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(size=14, color="#000000"), 
        legend.position="none",
        panel.background=element_blank(),
        plot.margin=unit(c(0.031,0,0.025,-0.2),"npc"))

# Create column of colored circles to represent phylogenetic clade
cladeLabels <- ggplot(data=denCladeData) +
  geom_point(mapping=aes(x=x, y=y, fill=cladeList), shape=21, color="black", stroke=1, size=4.5) +
  scale_fill_manual(values=c("#09c736", "#0dbecb", "#1127cf", "#9a16d3", "#d71a73", "#db601e", "#ffffff"), drop=FALSE) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=unit(c(0.001,0,0.087,-0.2),"npc"))

# Comconte dendrogram and heatmap into 1 figure
finalFigure <- plot_grid(contDenPlot, heatmap, cladeLabels, nrow=1, rel_widths=c(0.575, 0.4, 0.025))
print(finalFigure)