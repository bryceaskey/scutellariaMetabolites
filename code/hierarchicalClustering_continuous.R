library(tidyverse)
library(ggdendro)
library(cowplot)
library(viridis)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]
herbarium1_30 <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/preprocessed/20200214_herbarium1_30.csv")[, 2:6]
herbarium31_78 <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/preprocessed/20200812_herbarium31_78.csv")[, 2:6]
wrightii <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/preprocessed/20201007_wrightii.csv")[, 2:6]
cladeData <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/phylo-tree-clades.csv")

# Specify any species, organs, or metabolites to be removed
speciesToRemove <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastifolia", "hastafolia"), collapse = '|')
organsToRemove <- paste(c("flowers", "roots"), collapse = '|')
#metabolitesToRemove <- paste(c("chrysinG", "oroxyloside", "baicalin", "wogonoside", "acetoside", "apigeninG", "scutellarin", "hispidulinG"), collapse = '|')

# Remove barbata from fresh data - use only KR data
fresh <- fresh %>%
  filter(!grepl("barbata", species))

# Processing for fresh samples ----
freshData <- rbind(fresh, frozenKR, wrightii)

freshData$species <- factor(freshData$species)
freshData$organ <- factor(freshData$organ)
freshData$metabolite <- factor(freshData$metabolite)

freshData <- freshData %>%
  filter(!grepl(speciesToRemove, species)) %>%
  filter(!grepl(organsToRemove, organ)) #%>%
#filter(!grepl(metaboliteToRemove, metabolite))

# Fix naming errors
freshData$species <- as.character(freshData$species)
freshData$species[freshData$species=="RNA Seq"] <- "racemosa"
freshData$species[freshData$species=="havenesis"] <- "havanensis"
freshData$species[freshData$species=="hastafolia"] <- "hastifolia"
freshData$species[freshData$species=="pekinesis"] <- "pekinensis var. alpina"
freshData$species[freshData$species=="siphocampuloides"] <- "siphocampyloides"
freshData$species[freshData$species=="indica"] <- "indica var. coccinea"
freshData$species[freshData$species=="angustifolia ssp. angustifolia"] <- "angustifolia"
freshData$species[freshData$species=="drumondii"] <- "drummondii"
freshData$species[freshData$species=="holmgrenierum"] <- "holmgreniorum"
freshData$species[freshData$species=="leptosiplonsipkon"] <- "leptosiphon"
freshData$species[freshData$species=="multicularis"] <- "multicaulis"
freshData$species[freshData$species=="suffrutscens"] <- "suffrutescens"
freshData$species <- factor(freshData$species)

# Average together duplicate species, and leaf and shoot data
freshData <- freshData %>%
  group_by(species, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm)) %>%
  ungroup()

# Adjust fresh ppm to correct for dilution
# Data is saved at 5000 ppm. Divide by 5 to calculate at 1000 ppm (= umol/1 g FW)
freshData <- freshData %>%
  transmute(
    species=species,
    metabolite=metabolite,
    concentration_ppm=concentration_ppm/5,
    stError_ppm=stError_ppm/5
  )

# Processing for herbarium samples ----
herbariumData <- rbind(herbarium1_30, herbarium31_78)
herbariumData$species <- factor(herbariumData$species)
herbariumData$organ <- factor(herbariumData$organ)
herbariumData$metabolite <- factor(herbariumData$metabolite)

herbariumData <- herbariumData %>%
  filter(!grepl(speciesToRemove, species)) %>%
  filter(!grepl(organsToRemove, organ)) #%>%
#filter(!grepl(metaboliteToRemove, metabolite))

# Fix naming errors
herbariumData$species <- as.character(herbariumData$species)
herbariumData$species[herbariumData$species=="RNA Seq"] <- "racemosa"
herbariumData$species[herbariumData$species=="havenesis"] <- "havanensis"
herbariumData$species[herbariumData$species=="hastafolia"] <- "hastifolia"
herbariumData$species[herbariumData$species=="pekinesis"] <- "pekinensis var. alpina"
herbariumData$species[herbariumData$species=="siphocampuloides"] <- "siphocampyloides"
herbariumData$species[herbariumData$species=="indica"] <- "indica var. coccinea"
herbariumData$species[herbariumData$species=="angustifolia ssp. angustifolia"] <- "angustifolia"
herbariumData$species[herbariumData$species=="drumondii"] <- "drummondii"
herbariumData$species[herbariumData$species=="holmgrenierum"] <- "holmgreniorum"
herbariumData$species[herbariumData$species=="leptosiplonsipkon"] <- "leptosiphon"
herbariumData$species[herbariumData$species=="multicularis"] <- "multicaulis"
herbariumData$species[herbariumData$species=="suffrutscens"] <- "suffrutescens"
herbariumData$species <- factor(herbariumData$species)

# Average together duplicate species, and leaf and shoot data
herbariumData <- herbariumData %>%
  group_by(species, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm)) %>%
  ungroup()

# Adjust herbarium ppm to correct for dilution
# Data is saved at 1000 ppm. Divide by 10 to calculate at 100 ppm (= umol/0.1 g DW)
herbariumData <- herbariumData %>%
  transmute(
    species=species,
    metabolite=metabolite,
    concentration_ppm=concentration_ppm/10,
    stError_ppm=stError_ppm/10
  )

# Merge fresh and herbarium data ----
# If both fresh and herbarium samples are available for a species, the herbarium data should be used
# Iterate through herbarium species, and delete any rows in freshData which match
for(species in levels(herbariumData$species)){
  freshData <- freshData[!freshData$species==species, ]
}

# Combine fresh and herbarium data into a single dataframe
freshData$species <- as.character(freshData$species)
herbariumData$species <- as.character(herbariumData$species)
heatmapData <- rbind(freshData, herbariumData)

# Define function to convert units of ppm to micromol/L ----
ppm2microM <- function(input_ppm, metaboliteName){
  if(!is.na(input_ppm)){
    if(metaboliteName=="acetoside"){ #PubChem CID: 5281800 
      output_microM <- (input_ppm/624.6)*1000
    }else if(metaboliteName=="apigenin"){ #PubChem CID: 5280443
      output_microM <- (input_ppm/270.24)*1000
    }else if(metaboliteName=="apigeninG"){ #PubChem CID: 5319484
      output_microM <- (input_ppm/446.4)*1000
    }else if(metaboliteName=="baicalein"){ #PubChem CID: 5281605
      output_microM <- (input_ppm/270.24)*1000
    }else if(metaboliteName=="baicalin"){ #PubChem CID: 64982
      output_microM <- (input_ppm/446.4)*1000
    }else if(metaboliteName=="chrysin"){ #PubChem CID: 5281607
      output_microM <- (input_ppm/254.24)*1000
    }else if(metaboliteName=="chrysinG"){ #PubChem CID: 44257628
      output_microM <- (input_ppm/430.4)*1000
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
#speciesData <- subset(heatmapData, select=-c(concentration_ppm, stError_ppm, stError_microM))
#speciesData <- speciesData %>%
#  pivot_wider(names_from=metabolite, values_from=concentration_microM) %>%
#  remove_rownames %>%
#  column_to_rownames(var="species")
#speciesData <- scale(speciesData)

# Calculate distance matrix and perform heirarchical clustering of species
#speciesDist <- dist(speciesData, method="euclidian")
#speciesCluster <- hclust(d=speciesDist, method="average")
#print(paste("Correlation between cophenetic distances and actual distances:", cor(speciesDist, cophenetic(speciesCluster))))
#speciesDenData <- dendro_data(as.dendrogram(speciesCluster), type="rectangle")

# Transpose distance matrix and perform heirarchical clustering of flavonoids
#flavonoidData <- t(speciesData)
#flavonoidDist <- dist(flavonoidData, method="euclidian")
#flavonoidCluster <- hclust(d=flavonoidDist, method="average")
#flavonoidDenData <- dendro_data(as.dendrogram(flavonoidCluster), type="rectangle")

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
heatmapData$metabolite <- as.character(heatmapData$metabolite)
heatmapData$metabolite[heatmapData$metabolite=="acetoside"] <- "acteoside"
heatmapData$metabolite[heatmapData$metabolite=="apigeninG"] <- "apigenin 7-G"
heatmapData$metabolite[heatmapData$metabolite=="chrysinG"] <- "chrysin 7-G"
heatmapData$metabolite[heatmapData$metabolite=="hispidulinG"] <- "hispiduloside"
heatmapData$metabolite[heatmapData$metabolite=="oroxylinA"] <- "oroxylin A"
heatmapData$metabolite <- factor(heatmapData$metabolite)
flavonoidOrder <- c("apigenin", "apigenin 7-G", "scutellarein", "scutellarin", "hispidulin", "hispiduloside",
                    "chrysin", "chrysin 7-G", "baicalein", "baicalin", "oroxylin A", "oroxyloside", "wogonin", "wogonoside", "acteoside")
heatmapData$metabolite <- factor(heatmapData$metabolite, levels=flavonoidOrder)

heatmapData$metabolite_group <- NA
heatmapData$metabolite_group[heatmapData$metabolite %in% c("apigenin", "apigenin 7-G", "scutellarein", "scutellarin", "hispidulin", "hispiduloside")] <- "4'-hydroxyflavones"
heatmapData$metabolite_group[heatmapData$metabolite %in% c("chrysin", "chrysin 7-G", "baicalein", "baicalin", "oroxylin A", "oroxyloside", "wogonin", "wogonoside")] <- "4'-deoxyflavones"
heatmapData$metabolite_group[heatmapData$metabolite=="acteoside"] <- " "
heatmapData$metabolite_group <- factor(heatmapData$metabolite_group, levels=c("4'-hydroxyflavones", "4'-deoxyflavones", " "))

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
heatmapCladeData <- data.frame(x=1, y=1:length(speciesList), speciesList, cladeList)
heatmapCladeData$cladeList <- factor(heatmapCladeData$cladeList, levels=c(1, 2, 3, 4, 5))

# Add "S." to beginning of each species name
#speciesDenData$labels$label <- paste("S.", speciesDenData$labels$label)
heatmapData$species <- as.character(heatmapData$species)
heatmapData$species <- paste("S.", heatmapData$species)
heatmapData$species <- factor(heatmapData$species, levels=paste("S.", speciesList))
heatmapData$species <- fct_rev(heatmapData$species)

heatmapCladeData$speciesList <- paste("S.", heatmapCladeData$speciesList)
heatmapCladeData$speciesList <- factor(heatmapCladeData$speciesList, levels=levels(heatmapData$species))
heatmapCladeData$y <- (nrow(heatmapCladeData)+1)-heatmapCladeData$y

# Create column of colored circles to represent phylogenetic clade
cladeLabels <- ggplot(data=heatmapCladeData) +
  geom_point(mapping=aes(x=x, y=y, fill=cladeList), shape=21, color="black", size=4.5) +
  scale_fill_manual(values=c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#9632B8"), drop=FALSE, na.value=NA) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=margin(21,0,64,-610,"pt"))

# Label fresh data with asterisks
freshLabelData <- data.frame(x=numeric(), y=numeric(), speciesList=factor())
for(species in paste("S.", unique(freshData$species))){
  freshLabelData <- rbind(freshLabelData, heatmapCladeData[heatmapCladeData$speciesList==species, 1:3])
}
freshLabels <- ggplot(data=freshLabelData) + 
  geom_point(mapping=aes(x=x, y=y), shape=8, color="black", size=2.5) +
  ylim(1, length(speciesList)) +
  theme_void() +
  theme(plot.margin=margin(21,0,64,-49,"pt"))

# Capitalize first letter of each flavonoid name
capString <- function(string) {
  c <- strsplit(string, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
}
#flavonoidDenData$labels$label <- sapply(flavonoidDenData$labels$label, capString)
heatmapData$metabolite <- as.character(heatmapData$metabolite)
heatmapData$metabolite <- sapply(heatmapData$metabolite, capString)
heatmapData$metabolite <- factor(heatmapData$metabolite, levels=sapply(flavonoidOrder, capString))

# Create dendrogram for species
#speciesDenPlot <- ggplot() +
#  geom_segment(data=segment(speciesDenData), mapping=aes(x=x, y=y, xend=xend, yend=yend)) +
#  geom_text(data=label(speciesDenData), mapping=aes(x=x, y=y, label=label, hjust=1), nudge_y=6.5, size=4, fontface="italic") +
#  coord_flip() +
#  scale_y_reverse(expand=c(0.3, 0)) +
#  theme_dendro()

# Create dendrogram for flavonoids
#flavonoidDenPlot <- ggplot() +
#  geom_segment(data=segment(flavonoidDenData), mapping=aes(x=x, y=y, xend=xend, yend=yend)) +
#  geom_text(data=label(flavonoidDenData), mapping=aes(x=x, y=y, label=label), hjust=1, nudge_y=10, angle=90, size=2.25) +
#  theme_dendro() +
#  scale_y_reverse()

# Create heatmap
heatmap <- ggplot(data=heatmapData) +
  geom_raster(mapping=aes(x=species, y=metabolite, fill=concentration_microM)) +
  scale_fill_viridis() +
  labs(y="Flavonoid", fill="Concentration \n(µmol/g FW)") +
  coord_flip() +
  facet_grid(~metabolite_group, scales="free_x", space="free_x", switch="x") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12, margin=margin(5,0,0,0), color="black"),
        axis.title.y=element_blank(), axis.text.y=element_text(size=12, margin=margin(0,20,0,0), face="italic", color="black"),
        legend.position="top", legend.direction="horizontal", legend.title=element_text(size=14), legend.text=element_text(size=12), legend.key.size=unit(0.75, "cm"),
        panel.background=element_blank(), panel.spacing=unit(0, "lines"),
        strip.text.x=element_text(size=12, color="black"), strip.placement="outside", strip.background=element_rect(fill="white"), axis.title=element_blank())

# Combine dendrograms and heatmap into 1 figure
heatmap <- plot_grid(heatmap, cladeLabels, freshLabels, nrow=1, rel_widths=c(1.5, 0.05, 0.05))
#speciesDendrogram <- plot_grid(speciesDenPlot, cladeLabels, nrow=1)
#flavonoidDendogram <- plot_grid(flavonoidDenPlot)

# Export dendrograms and heatmaps separately
ggsave(filename="C:/Users/Bryce/Research/scutellariaMetabolites/figures/heatmaps/heatmap.png",
  plot=heatmap,
  device=png(),
  width=18, height=45, units="cm")
dev.off()

#ggsave(filename="C:/Users/Bryce/Documents/scutellariaMetabolites/figures/heatmaps/speciesDendrogram.png",
#  plot=speciesDendrogram,
#  device=png(),
#  width=30, height=30, units="cm")

#ggsave(filename="C:/Users/Bryce/Documents/scutellariaMetabolites/figures/heatmaps/flavonoidDendrogram.png",
#  plot=flavonoidDendogram,
#  device=png(),
#  width=10, height=4, units="cm")

# Print summary statistics
print(paste("# of species included:", nrow(heatmapCladeData)))
print(paste("# of species assigned to a clade:", sum(!is.na(heatmapCladeData$cladeList))))
flavonoidAbundance <- data.frame(flavonoid=character(), totalCount=numeric(), totalPct=numeric(), avgConc=numeric(),
                                 clade1Pct=numeric(), clade2Pct=numeric(), clade3Pct=numeric(), clade4Pct=numeric(), clade5Pct=numeric())
for(clade in 1:5){
  assign(paste("clade", clade, "Sp", sep=""), as.character(heatmapCladeData$speciesList[heatmapCladeData$cladeList==clade & !is.na(heatmapCladeData$cladeList)]))
  assign(paste("clade", clade, "heatmapData", sep=""), heatmapData[heatmapData$species %in% get(paste("clade", clade, "Sp", sep="")), ])
}
for(flavonoid in levels(heatmapData$metabolite)){
  occurrence <- data.frame(flavonoid=flavonoid, 
                           totalCount=sum(heatmapData$metabolite==flavonoid & heatmapData$concentration_microM>0),
                           totalPct=sum(heatmapData$metabolite==flavonoid & heatmapData$concentration_microM>0)/nrow(heatmapCladeData),
                           avgConc=mean(heatmapData$concentration_microM[heatmapData$metabolite==flavonoid & heatmapData$concentration_microM>0]),
                           clade1Pct=sum(clade1heatmapData$metabolite==flavonoid & clade1heatmapData$concentration_microM>0)/sum(clade1heatmapData$metabolite==flavonoid),
                           clade2Pct=sum(clade2heatmapData$metabolite==flavonoid & clade2heatmapData$concentration_microM>0)/sum(clade2heatmapData$metabolite==flavonoid),
                           clade3Pct=sum(clade3heatmapData$metabolite==flavonoid & clade3heatmapData$concentration_microM>0)/sum(clade3heatmapData$metabolite==flavonoid),
                           clade4Pct=sum(clade4heatmapData$metabolite==flavonoid & clade4heatmapData$concentration_microM>0)/sum(clade4heatmapData$metabolite==flavonoid),
                           clade5Pct=sum(clade5heatmapData$metabolite==flavonoid & clade5heatmapData$concentration_microM>0)/sum(clade5heatmapData$metabolite==flavonoid))
  flavonoidAbundance <- rbind(flavonoidAbundance, occurrence)
}
flavonoidAbundance <- flavonoidAbundance[order(flavonoidAbundance$totalCount), ]
print("Flavonoids in order of abundance:")
print(flavonoidAbundance)

deoxyflavones <- heatmapData[heatmapData$metabolite %in% c("Chrysin", "Chrysin 7-G", "Baicalin", "Baicalein", "Wogonin", "Wogonoside", "Oroxylin A", "Oroxyloside"), ]
deoxyflavones <- deoxyflavones %>%
  group_by(species) %>%
  summarise(totalDeoxyflavoneContent=sum(concentration_microM))
print(paste("Species with no deoxyflavones detected:", sum(deoxyflavones$totalDeoxyflavoneContent==0)))