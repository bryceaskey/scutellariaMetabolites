library(tidyverse)
library(ggdendro)
library(cowplot)
library(viridis)
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(mixtools)
library(robustbase)
library(ggrepel)
library(ggsci)

# Generate heatmap ----

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20190813_fresh.csv")
frozenKR <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200117_frozenKR.csv")
herbarium1_30 <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200214_herbarium1_30.csv")
herbarium31_78 <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200812_herbarium31_78.csv")
wrightii <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201007_wrightii.csv")
suffrutescencs <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201119_suffrutescens.csv")
cladeData <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/herbarium/phylo-tree-clades.csv")

# Remove barbata from frozenKR data - use only fresh data
frozenKR <- frozenKR %>%
  filter(!grepl("barbata", species))

# Processing for fresh samples
freshData <- rbind(fresh, frozenKR, wrightii)
freshData$species <- factor(freshData$species)
freshData$organ <- factor(freshData$organ)
freshData$metabolite <- factor(freshData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa_071119", "racemosa_MS", "racemosa_SC", "hastifolia", "arenicola", "havanensis"), collapse = '|')
excludeOrgans <- paste(c("flowers", "roots"), collapse = '|')
excludeMetabolites <- paste(c("isoscutellarin"), collapse = '|')
freshData <- freshData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) %>%
  filter(!grepl(excludeMetabolites, metabolite))

# Fix naming errors
freshData$species <- as.character(freshData$species)
freshData$species[freshData$species=="racemosa_RNAseq"] <- "racemosa"
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

# Processing for herbarium samples
herbariumData <- rbind(herbarium1_30, herbarium31_78)
herbariumData$species <- factor(herbariumData$species)
herbariumData$organ <- factor(herbariumData$organ)
herbariumData$metabolite <- factor(herbariumData$metabolite)

herbariumData <- herbariumData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) %>%
  filter(!grepl(excludeMetabolites, metabolite))

# Fix naming errors
herbariumData$species <- as.character(herbariumData$species)
herbariumData$species[herbariumData$species=="racemosa_RNAseq"] <- "racemosa"
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

# Merge fresh and herbarium data
# If both fresh and herbarium samples are available for a species, the herbarium data should be used
# Iterate through herbarium species, and delete any rows in freshData which match
for(species in levels(herbariumData$species)){
  freshData <- freshData[!freshData$species==species, ]
}

# Combine fresh and herbarium data into a single dataframe
freshData$species <- as.character(freshData$species)
herbariumData$species <- as.character(herbariumData$species)
heatmapData <- rbind(freshData, herbariumData)

# Define function to convert units of ppm to micromol/L
ppm2microM <- function(input_ppm, metaboliteName){
  if(!is.na(input_ppm)){
    if(metaboliteName=="acteoside"){ #PubChem CID: 5281800 
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
    }else if(metaboliteName=="isoscutellarin"){
      output_microM <- (input_ppm/462.4)*1000
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

# Adjust species order in heatmap data to match order in phylogenetic tree
speciesOrder <- cladeData$species
heatmapData$species <- factor(heatmapData$species, levels=speciesOrder)
heatmapData <- heatmapData[order(heatmapData$species), ]
heatmapData$species <- droplevels(heatmapData$species)

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

# Assign metabolites to groups to generated faceted heatmap
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
heatmapData$species <- as.character(heatmapData$species)
heatmapData$species <- paste("S.", heatmapData$species)
heatmapData$species <- factor(heatmapData$species, levels=paste("S.", speciesList))
heatmapData$species <- fct_rev(heatmapData$species)

heatmapCladeData$speciesList <- paste("S.", heatmapCladeData$speciesList)
heatmapCladeData$speciesList <- factor(heatmapCladeData$speciesList, levels=levels(heatmapData$species))
heatmapCladeData$y <- (nrow(heatmapCladeData)+1)-heatmapCladeData$y

# Create column of colored circles to represent phylogenetic clade
cladeLabels <- ggplot(data=heatmapCladeData) +
  geom_point(mapping=aes(x=x, y=y, fill=cladeList), shape=21, color="black", size=2.5) +
  scale_fill_manual(values=c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#9632B8"), drop=FALSE, na.value=NA) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=margin(4.75,0,37,-288,"pt"))

# Label fresh data with asterisks
freshLabelData <- data.frame(x=numeric(), y=numeric(), speciesList=factor())
for(species in paste("S.", unique(freshData$species))){
  freshLabelData <- rbind(freshLabelData, heatmapCladeData[heatmapCladeData$speciesList==species, 1:3])
}
freshLabels <- ggplot(data=freshLabelData) + 
  geom_point(mapping=aes(x=x, y=y), shape=8, color="black", size=1) +
  ylim(1, length(speciesList)) +
  theme_void() +
  theme(plot.margin=margin(4.75,0,37,-30,"pt"))

# Capitalize first letter of each flavonoid name
capString <- function(string) {
  c <- strsplit(string, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
}
#flavonoidDenData$labels$label <- sapply(flavonoidDenData$labels$label, capString)
heatmapData$metabolite <- as.character(heatmapData$metabolite)
heatmapData$metabolite <- sapply(heatmapData$metabolite, capString)
heatmapData$metabolite <- factor(heatmapData$metabolite, levels=sapply(flavonoidOrder, capString))

# Create heatmap
heatmap <- ggplot(data=heatmapData) +
  geom_raster(mapping=aes(x=species, y=metabolite, fill=concentration_microM)) +
  scale_fill_viridis() +
  labs(y="Flavonoid", fill="Concentration \n(µmol/g FW)") +
  coord_flip() +
  facet_grid(~metabolite_group, scales="free_x", space="free_x", switch="x") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=6, margin=margin(2,0,-5,0), color="black"),
        axis.title.y=element_blank(), axis.text.y=element_text(size=6, margin=margin(0,10,0,0), face="italic", color="black"),
        legend.position="top", legend.direction="horizontal", legend.title=element_text(size=8),
        legend.text=element_text(size=6), legend.key.size=unit(0.3, "cm"), legend.margin=margin(0,0,0,0), legend.box.margin=margin(b=-8),
        panel.background=element_blank(), panel.spacing=unit(0, "lines"),
        strip.text.x=element_text(size=6, color="black"), strip.placement="outside", strip.background=element_rect(fill="white"), axis.title=element_blank(),
        axis.ticks=element_line(color="black"))

# Combine heatmap, clade labels, and fresh data labels into 1 figure
heatmap <- plot_grid(heatmap, cladeLabels, freshLabels, nrow=1, rel_widths=c(1.5, 0.05, 0.05))
print(heatmap)

# Generate MCA plots ----
# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20190813_fresh.csv")
frozenKR <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200117_frozenKR.csv")
herbarium1_30 <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200214_herbarium1_30.csv")
herbarium31_78 <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20200812_herbarium31_78.csv")
wrightii <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201007_wrightii.csv")
cladeData <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/herbarium/phylo-tree-clades.csv")

# Processing for fresh samples
freshData <- rbind(fresh, frozenKR, wrightii)

freshData$species <- factor(freshData$species)
freshData$organ <- factor(freshData$organ)
freshData$metabolite <- factor(freshData$metabolite)

freshData <- freshData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeMetabolites, organ)) %>%
  filter(!grepl(excludeMetabolites, metabolite))

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

# Processing for herbarium samples
herbariumData <- rbind(herbarium1_30, herbarium31_78)
herbariumData$species <- factor(herbariumData$species)
herbariumData$organ <- factor(herbariumData$organ)
herbariumData$metabolite <- factor(herbariumData$metabolite)

herbariumData <- herbariumData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) %>%
  filter(!grepl(excludeMetabolites, metabolite))

print("error here")

# Fix naming errors
herbariumData$species <- as.character(herbariumData$species)
herbariumData$species[herbariumData$species=="racemosa_RNASeq"] <- "racemosa"
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

# Merge fresh and herbarium data
# If both fresh and herbarium samples are available for a species, the herbarium data should be used
# Iterate through herbarium species, and delete any rows in freshData which match
for(species in levels(herbariumData$species)){
  freshData <- freshData[!freshData$species==species, ]
}

# Combine fresh and herbarium data into a single dataframe
freshData$species <- as.character(freshData$species)
herbariumData$species <- as.character(herbariumData$species)
allData <- rbind(freshData, herbariumData)
allData$species <- factor(allData$species)

# Function to convert units of ppm to micromol/L
ppm2microM <- function(input_ppm, metaboliteName){
  if(!is.na(input_ppm)){
    if(metaboliteName=="acteoside"){ #PubChem CID: 5281800 
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
    }else if(metaboliteName=="isoscutellarin"){
      output_microM <- (input_ppm/462.4)*1000
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
concentration_microM <- vector(mode="numeric", length=nrow(allData))
stError_microM <- vector(mode="numeric", length=nrow(allData))
for(i in 1:nrow(allData)){
  concentration_microM[i] <- ppm2microM(allData$concentration_ppm[i], allData$metabolite[i])
  stError_microM[i] <- ppm2microM(allData$stError_ppm[i], allData$metabolite[i])
}
allData$concentration_microM <- concentration_microM
allData$stError_microM <- stError_microM

allData$metabolite <- as.character(allData$metabolite)
allData$metabolite[allData$metabolite=="acetoside"] <- "acteoside"
allData$metabolite[allData$metabolite=="acetoside"] <- "acteoside"
allData$metabolite[allData$metabolite=="apigeninG"] <- "apigenin 7-G"
allData$metabolite[allData$metabolite=="chrysinG"] <- "chrysin 7-G"
allData$metabolite[allData$metabolite=="hispidulinG"] <- "hispiduloside"
allData$metabolite[allData$metabolite=="oroxylinA"] <- "oroxylin A"
allData$metabolite <- factor(allData$metabolite)

# Transform data into wide format to use for heirarchical clustering 
speciesData <- subset(allData, select=-c(concentration_ppm, stError_ppm, stError_microM))
speciesData <- speciesData %>%
  pivot_wider(names_from=metabolite, values_from=concentration_microM) %>%
  remove_rownames %>%
  column_to_rownames(var="species")

# Capitalize first letter of flavonoid names
capString <- function(string) {
  c <- strsplit(string, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
}
colnames(speciesData) <- sapply(colnames(speciesData), capString)

# Add column to indicate clade which the species belongs to
cladeList <- vector(mode="numeric", length=nrow(speciesData))
for (i in 1:nrow(speciesData)){
  species <- rownames(speciesData)[i]
  clade <- cladeData$clade[cladeData$species==species]
  if (length(clade)>0 && !is.na(clade)){
    cladeList[i] <- clade
  }else{
    cladeList[i] <- NA
  }
} 
speciesData$clade <- factor(cladeList)

# For MCA with binary data
for(i in 1:15){
  speciesData[, i] <- as.logical(speciesData[, i])
}
speciesData <- speciesData[c("Apigenin", "Apigenin 7-G", "Scutellarein", "Scutellarin", "Hispidulin", "Hispiduloside",
                             "Chrysin", "Chrysin 7-G", "Baicalein", "Baicalin", "Oroxylin A", "Oroxyloside", "Wogonin", "Wogonoside", "Acteoside",
                             "clade")]
pca_data <- MCA(speciesData[, c(1:15)], graph=FALSE)

# Extract data for plotting
pca_inds <- data.frame(pc1_ind=pca_data$ind$coord[,1],
                       pc2_ind=pca_data$ind$coord[,2],
                       clade=speciesData$clade)
pc1_expl <- round(pca_data$eig[1,2],2)
pc2_expl <- round(pca_data$eig[2,2],2)

# Create 95% and 80% confidence ellipses
ellipse_allClade_0.95 <- data.frame(pc1=numeric(), pc2=numeric(), clade=factor())
ellipse_allClade_0.80 <- data.frame(pc1=numeric(), pc2=numeric(), clade=factor())
for(i in levels(pca_inds$clade)){
  covariance <- covMcd(filter(pca_inds, clade==i)[,1:2], alpha=0.95)
  ellipse_center <- unname(covariance$center)
  ellipse_cov <- covariance$cov
  
  ellipse_matrix_0.95 <- ellipse(ellipse_center, ellipse_cov, draw=FALSE, alpha=0.05)
  ellipse_df_0.95 <- data.frame(pc1=ellipse_matrix_0.95[,1], pc2=ellipse_matrix_0.95[,2], clade=i)
  ellipse_allClade_0.95 <- rbind(ellipse_allClade_0.95, ellipse_df_0.95)
  
  ellipse_matrix_0.80 <- ellipse(ellipse_center, ellipse_cov, draw=FALSE, alpha=0.20)
  ellipse_df_0.80 <- data.frame(pc1=ellipse_matrix_0.80[,1], pc2=ellipse_matrix_0.80[,2], clade=i)
  ellipse_allClade_0.80 <- rbind(ellipse_allClade_0.80, ellipse_df_0.80)
}

pcaPlot <- ggplot() +
  geom_polygon(data=ellipse_allClade_0.80, mapping=aes(x=pc1, y=pc2, color=clade, group=clade), linetype=1, fill=NA, size=1) +
  geom_hline(mapping=aes(yintercept=0), color="darkgray", linetype="dashed") +
  geom_vline(mapping=aes(xintercept=0), color="darkgray", linetype="dashed") +
  geom_point(data=pca_inds, mapping=aes(x=pc1_ind, y=pc2_ind, fill=clade), shape=21, color="black", size=2.5, position=position_jitter(h=0.05, w=0.05)) +
  scale_fill_manual(values=c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#9632B8"), drop=FALSE, na.value=NA) +
  scale_color_manual(values=c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#9632B8"), drop=FALSE, na.value=NA) +
  coord_fixed(ratio=1) +
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1, 1)) +
  xlab(paste("PC1 (", pc1_expl, "%)", sep="")) +
  ylab(paste("PC2 (", pc2_expl, "%)", sep="")) +
  theme_classic() +
  theme(legend.position="none",
        panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
        axis.title=element_text(size=8), axis.text=element_text(size=6),
        text=element_text(color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        plot.margin=margin(t=20),
        axis.ticks=element_line(color="black")
        )
print(pcaPlot)

pca_vars <- data.frame(get_mca_var(pca_data, "var")$coord[,1:2])
pca_vars <- rownames_to_column(pca_vars, var="variable")
pca_vars$metaboliteClass <- c(rep("4'-hydroxyflavone", 12), rep("4'-deoxyflavone", 16), rep("Acteoside", 2))
pca_vars$metaboliteClass <- factor(pca_vars$metaboliteClass, levels=c("4'-hydroxyflavone", "4'-deoxyflavone", "Acteoside"))
varPlot <- ggplot(data=pca_vars, mapping=aes(x=Dim.1, y=Dim.2, label=variable, fill=metaboliteClass)) +
  geom_hline(mapping=aes(yintercept=0), color="darkgray", linetype="dashed") +
  geom_vline(mapping=aes(xintercept=0), color="darkgray", linetype="dashed") +
  geom_point(size=2.5, shape=21, color="black") +
  geom_text_repel(data=subset(pca_vars, Dim.1<0), color="black", nudge_x=-0.5, direction="y", point.padding=0.6, box.padding=0.15, min.segment.length=0.1, size=6/.pt) +
  geom_text_repel(data=subset(pca_vars, Dim.1>0), color="black", nudge_x=0.5, direction="y", point.padding=0.6, box.padding=0.15, min.segment.length=0.1, size=6/.pt) +
  scale_fill_npg(name="")+
  coord_fixed(ratio=1) +
  coord_cartesian(xlim=c(-1.55, 1.55), ylim=c(-1.55, 1.55)) +
  xlab(paste("PC1 (", pc1_expl, "%)", sep="")) +
  ylab(paste("PC2 (", pc2_expl, "%)", sep="")) +
  scale_x_continuous(breaks=c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5)) +
  scale_y_continuous(breaks=c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5)) +
  theme_classic() +
  theme(legend.position="top",
        panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
        axis.title=element_text(size=8), axis.text=element_text(size=6),
        legend.title=element_text(size=8), legend.text=element_text(size=8, margin=margin(l=-5, r=5)), legend.margin=margin(0,0,0,0), legend.box.margin=margin(t=40),
        text=element_text(color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        plot.margin=margin(b=20),
        axis.ticks=element_line(color="black")
        )
print(varPlot)

pcaVarPlot <- plot_grid(pcaPlot, varPlot, ncol=1, nrow=2, rel_heights=c(1,1.3), labels=c("B", "C"), hjust=0.5, vjust=c(1.5, 4.5))

# Merge plots and export as pdf ----
totalPlot <- plot_grid(heatmap, pcaVarPlot, ncol=2, nrow=1, rel_widths=c(1,1.1), labels=c("A"))


# Export figure
ggsave(filename="C:/Users/Bryce/Research/scutellariaMetabolites/figures/heatmapPCA/20210724_noIsoscutellarin.pdf",
      plot=totalPlot,
      device=pdf(),
      width=7.25, height=8.5, units="in")
dev.off()

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