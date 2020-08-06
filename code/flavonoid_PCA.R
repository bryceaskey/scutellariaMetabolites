library(tidyverse)
library(ggplot2)
library(grid)
library(dplyr)
library(FactoMineR)

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

for(i in 1:15){
  speciesData[, i] <- as.logical(speciesData[, i])
}

pca_data <- MCA(speciesData[, c(1:15)], ncp=2, graph=TRUE)
pc1_ind <- pca_data$ind$coord[,1]
pc2_ind <- pca_data$ind$coord[,2]
pc1_expl <- round(pca_data$eig[1,2],2)
pc2_expl <- round(pca_data$eig[2,2],2)
clade <- speciesData$clade

pcaPlot <- ggplot() +
  geom_point(aes(x=pc1_ind, y=pc2_ind, fill=clade), color="black", pch=21, size=7) +
  scale_fill_manual(values=c("#62e8ec", "#90dfb0", "#c6ce86", "#f0b682", "#ffa2a2", "#FFFFFF")) +
  #coord_fixed(ratio=1) +
  #xlab(paste("PC1 (", pc1_expl, "%)", sep="")) +
  #ylab(paste("PC2 (", pc2_expl, "%)", sep="")) +
  xlab("Principal component 1") +
  ylab("Principal component 2") +
  theme_classic() +
  theme(legend.position="none",
        panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
        axis.title=element_text(size=14), axis.text=element_text(size=12))
print(pcaPlot)