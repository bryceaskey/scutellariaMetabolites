library(tidyverse)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(viridis)
library(dplyr)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]
herbarium1_30 <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200214_herbarium1_30.csv")[, 2:6]
cladeData <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/phylo-tree-clades.csv")

# Adjust herbarium ppm to correct for dilution
herbarium1_30 <- herbarium1_30 %>%
  transmute(
    species=species,
    organ=organ,
    metabolite=metabolite,
    concentration_ppm=concentration_ppm/2,
    stError_ppm=stError_ppm/2
  )

# Combine all data into a single data frame and change classifiers (species, organs, metabolites)
# into factors
allData <- rbind(fresh, frozenKR, herbarium1_30)
allData$species <- as.factor(allData$species)
allData$organ <- as.factor(allData$organ)
allData$metabolite <- as.factor(allData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastifolia", "hastafolia"), collapse = '|')
excludeOrgans <- paste(c("flowers", "roots"), collapse = '|')
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
allData$species[allData$species=="pekinesis"] <- "pekinensis var. alpina"
allData$species[allData$species=="siphocampuloides"] <- "siphocampyloides"
allData$species[allData$species=="indica"] <- "indica var. coccinea"
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

# Transform data into wide format to use for heirarchical clustering 
speciesData <- subset(heatmapData, select=-c(concentration_ppm, stError_ppm, stError_microM))
speciesData <- speciesData %>%
  pivot_wider(names_from=metabolite, values_from=concentration_microM) %>%
  remove_rownames %>%
  column_to_rownames(var="species")
speciesData <- scale(speciesData)

# Calculate distance matrix and perform heirarchical clustering of species
speciesDist <- dist(speciesData, method="euclidian")
speciesCluster <- hclust(d=speciesDist, method="average")
print(paste("Correlation between cophenetic distances and actual distances:", cor(speciesDist, cophenetic(speciesCluster))))
speciesDenData <- dendro_data(as.dendrogram(speciesCluster), type="rectangle")

# Transpose distance matrix and perform heirarchical clustering of flavonoids
flavonoidData <- t(speciesData)
flavonoidDist <- dist(flavonoidData, method="euclidian")
flavonoidCluster <- hclust(d=flavonoidDist, method="average")
flavonoidDenData <- dendro_data(as.dendrogram(flavonoidCluster), type="rectangle")

# Adjust species order in heatmap data to match order in dendrogram
#speciesOrder <- label(speciesDenData)$label
# Adjust species order in heatmap data to match order in phylogenetic tree
speciesOrder <- cladeData$species
heatmapData$species <- factor(heatmapData$species, levels=speciesOrder)
heatmapData <- heatmapData[order(heatmapData$species), ]
heatmapData$species <- droplevels(heatmapData$species)

# Adjust flavonoid order in heatmap data to match order in dendrogram
#flavonoidOrder <- label(flavonoidDenData)$label
# Adjust flavonoid order in heatmap data to match order in biosynthetic pathway
flavonoidOrder <- c("apigenin", "apigeninG", "scutellarein", "scutellarin", "hispidulin", "hispidulinG",
                    "chrysin", "chrysinG", "baicalein", "baicalin", "oroxylinA", "oroxyloside", "wogonin", "wogonoside", "acetoside")
heatmapData$metabolite <- factor(heatmapData$metabolite, levels=flavonoidOrder)

# Load clade data from phylogenetic tree generated from chloroplast genome
speciesList <- vector(mode="character", length=length(levels(heatmapData$species)))
cladeList <- vector(mode="numeric", length=length(levels(heatmapData$species)))
for (i in 1:length(levels(heatmapData$species))){
  species <- levels(heatmapData$species)[i]
  speciesList[i] <- species
  clade <- cladeData$clade[cladeData$species==species]
  if (length(clade)>0 && !is.na(clade)){
    cladeList[i] <- clade
  }else{
    cladeList[i] <- NA
  }
} 
heatmapCladeData <- data.frame(x=1, y=1:nrow(label(speciesDenData)), speciesList, cladeList)
heatmapCladeData$cladeList <- factor(heatmapCladeData$cladeList, levels=c(1, 2, 3, 4, 5))

# Add "S." to beginning of each species name
speciesDenData$labels$label <- paste("S.", speciesDenData$labels$label)
heatmapData$species <- as.character(heatmapData$species)
heatmapData$species <- paste("S.", heatmapData$species)
heatmapData$species <- factor(heatmapData$species, levels=paste("S.", speciesList))
heatmapData$species <- fct_rev(heatmapData$species)

heatmapCladeData$speciesList <- paste("S.", heatmapCladeData$speciesList)
heatmapCladeData$speciesList <- factor(heatmapCladeData$speciesList, levels=levels(heatmapData$species))
heatmapCladeData$y <- 40-heatmapCladeData$y

# Create column of colored circles to represent phylogenetic clade
cladeLabels <- ggplot(data=heatmapCladeData) +
  geom_point(mapping=aes(x=x, y=y, fill=cladeList), shape=21, color="black", stroke=0.6, size=5.5) +
  scale_fill_manual(values=c("#62e8ec", "#90dfb0", "#c6ce86", "#f0b682", "#ffa2a2", "#FFFFFF"), drop=FALSE) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=margin(44,0,58,-715,"pt"))

# Capitalize first letter of each flavonoid name
capString <- function(string) {
  c <- strsplit(string, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
}
flavonoidDenData$labels$label <- sapply(flavonoidDenData$labels$label, capString)
heatmapData$metabolite <- as.character(heatmapData$metabolite)
heatmapData$metabolite <- sapply(heatmapData$metabolite, capString)
heatmapData$metabolite <- factor(heatmapData$metabolite, levels=sapply(flavonoidOrder, capString))

# Create dendrogram for species
speciesDenPlot <- ggplot() +
  geom_segment(data=segment(speciesDenData), mapping=aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(speciesDenData), mapping=aes(x=x, y=y, label=label, hjust=1), nudge_y=6.5, size=4, fontface="italic") +
  coord_flip() +
  scale_y_reverse(expand=c(0.3, 0)) +
  theme_dendro()

# Create dendrogram for flavonoids
flavonoidDenPlot <- ggplot() +
  geom_segment(data=segment(flavonoidDenData), mapping=aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(flavonoidDenData), mapping=aes(x=x, y=y, label=label), hjust=1, nudge_y=10, angle=90, size=2.25) +
  theme_dendro() +
  scale_y_reverse()

# Create heatmap without species or flavonoid labels
heatmap <- ggplot(data=heatmapData) +
  geom_raster(mapping=aes(x=species, y=metabolite, fill=concentration_microM)) +
  scale_fill_viridis() +
  labs(y="Flavonoid", fill=expression(paste("Conc (", mu, "M)", sep=""))) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12, margin=margin(5,0,0,0), color="black"),
        axis.title.y=element_blank(), axis.text.y=element_text(size=12, margin=margin(0,20,0,0), face="italic"),
        legend.position="top", legend.direction="horizontal", legend.title=element_text(size=14), legend.text=element_text(size=12), legend.key.size=unit(0.75, "cm"),
        panel.background=element_blank())

# Combine dendrograms and heatmap into 1 figure
heatmap <- plot_grid(heatmap, cladeLabels, nrow=1, rel_widths=c(1.5, 0.05))
#speciesDendrogram <- plot_grid(speciesDenPlot, cladeLabels, nrow=1)
flavonoidDendogram <- plot_grid(flavonoidDenPlot)

# Export dendrograms and heatmaps separately
ggsave(filename="C:/Users/Bryce/Documents/scutellariaMetabolites/figures/heatmaps/heatmap.png",
  plot=heatmap,
  device=png(),
  width=18, height=30, units="cm")

#ggsave(filename="C:/Users/Bryce/Documents/scutellariaMetabolites/figures/heatmaps/speciesDendrogram.png",
#  plot=speciesDendrogram,
#  device=png(),
#  width=30, height=30, units="cm")

ggsave(filename="C:/Users/Bryce/Documents/scutellariaMetabolites/figures/heatmaps/flavonoidDendrogram.png",
  plot=flavonoidDendogram,
  device=png(),
  width=10, height=4, units="cm")