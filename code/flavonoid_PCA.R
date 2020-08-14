library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpubr)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]
herbarium1_30 <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200214_herbarium1_30.csv")[, 2:6]
herbarium31_78 <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200812_herbarium31_78.csv")[, 2:6]
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

herbarium31_78 <- herbarium31_78 %>%
  transmute(
    species=species,
    organ=organ,
    metabolite=metabolite,
    concentration_ppm=concentration_ppm/2,
    stError_ppm=stError_ppm/2
  )

# Specify any species, organs, or metabolites to be removed
speciesToRemove <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastifolia", "hastafolia"), collapse = '|')
organsToRemove <- paste(c("flowers", "roots"), collapse = '|')
#metabolitesToRemove <- paste(c("chrysinG", "oroxyloside", "baicalin", "wogonoside", "acetoside", "apigeninG", "scutellarin", "hispidulinG"), collapse = '|')

# Processing for fresh samples ----
freshData <- rbind(fresh, frozenKR)

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
freshData$species <- factor(freshData$species)

# Average together duplicate species, and leaf and shoot data
freshData <- freshData %>%
  group_by(species, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm))

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
herbariumData$species <- factor(herbariumData$species)

# Average together duplicate species, and leaf and shoot data
herbariumData <- herbariumData %>%
  group_by(species, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm))

# Merge fresh and herbarium data ----
# If both fresh and herbarium samples are available for a species, the herbarium data should be used
# Iterate through herbarium species, and delete any rows in freshData which match
for(herbariumSpecies in levels(herbariumData$species)){
  freshData <- freshData[!freshData$species==herbariumSpecies, ]
}

# Combine fresh and herbarium data into a single dataframe
freshData$species <- as.character(freshData$species)
herbariumData$species <- as.character(herbariumData$species)
allData <- rbind(freshData, herbariumData)


# Function to convert units of ppm to micromol/L
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
concentration_microM <- vector(mode="numeric", length=nrow(allData))
stError_microM <- vector(mode="numeric", length=nrow(allData))
for(i in 1:nrow(allData)){
  concentration_microM[i] <- ppm2microM(allData$concentration_ppm[i], allData$metabolite[i])
  stError_microM[i] <- ppm2microM(allData$stError_ppm[i], allData$metabolite[i])
}
allData$concentration_microM <- concentration_microM
allData$stError_microM <- stError_microM

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
pca_data <- MCA(speciesData, quali.sup=16,  graph=FALSE)

# For PCA with continuous data
#pca_data <- PCA(speciesData[, c(1:15)], ncp=2, scale.unit=FALSE, graph=TRUE)

# Extract data for plotting
pca_inds <- data.frame(pc1_ind=pca_data$ind$coord[,1],
                       pc2_ind=pca_data$ind$coord[,2],
                       clade=speciesData$clade)
pc1_expl <- round(pca_data$eig[1,2],2)
pc2_expl <- round(pca_data$eig[2,2],2)

# Create 95% confidence ellipses
#clade1_cor <- cor(filter(pca_inds, clade==1)[,1:2])
#clade1_ellipse <- as.data.frame(ellipse(clade1_cor*pca_data$eig[1,1], centre=colMeans(filter(pca_inds, clade==1)[,1:2]), level=0.95))

#clade2_cor <- cor(filter(pca_inds, clade==2)[,1:2])
#clade2_ellipse <- as.data.frame(ellipse(clade2_cor*pca_data$eig[1,1], centre=colMeans(filter(pca_inds, clade==2)[,1:2]), level=0.95))

#clade3_cor <- cor(filter(pca_inds, clade==3)[,1:2])
#clade3_ellipse <- as.data.frame(ellipse(clade3_cor*pca_data$eig[1,1], centre=colMeans(filter(pca_inds, clade==3)[,1:2]), level=0.95))

#clade4_cor <- cor(filter(pca_inds, clade==4)[,1:2])
#clade4_ellipse <- as.data.frame(ellipse(clade4_cor*pca_data$eig[1,1], centre=colMeans(filter(pca_inds, clade==4)[,1:2]), level=0.95))

#clade5_cor <- cor(filter(pca_inds, clade==5)[,1:2])
#clade5_ellipse <- as.data.frame(ellipse(clade5_cor*pca_data$eig[1,1], centre=colMeans(filter(pca_inds, clade==5)[,1:2]), level=0.95))

pcaPlot <- ggplot() +
  #geom_polygon(data=clade1_ellipse, mapping=aes(x=pc1_ind, y=pc2_ind), fill=NA, color="#62e8ec") +
  #geom_polygon(data=clade2_ellipse, mapping=aes(x=pc1_ind, y=pc2_ind), fill=NA, color="#90dfb0") +
  #geom_polygon(data=clade3_ellipse, mapping=aes(x=pc1_ind, y=pc2_ind), fill=NA, color="#c6ce86") +
  #geom_polygon(data=clade4_ellipse, mapping=aes(x=pc1_ind, y=pc2_ind), fill=NA, color="#f0b682") +
  #geom_polygon(data=clade5_ellipse, mapping=aes(x=pc1_ind, y=pc2_ind), fill=NA, color="#ffa2a2") +
  geom_point(data=pca_inds, mapping=aes(x=pc1_ind, y=pc2_ind, fill=clade), color="black", pch=21, size=7) +
  scale_fill_manual(values=c("#62e8ec", "#90dfb0", "#c6ce86", "#f0b682", "#ffa2a2", "#FFFFFF")) +
  coord_fixed(ratio=1) +
  xlab(paste("PC1 (", pc1_expl, "%)", sep="")) +
  ylab(paste("PC2 (", pc2_expl, "%)", sep="")) +
  theme_classic() +
  theme(legend.position="none",
        panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
        axis.title=element_text(size=14), axis.text=element_text(size=12))
print(pcaPlot)