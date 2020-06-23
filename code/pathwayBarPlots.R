# Code to create metabolite bar graphs for pathway figure

library(ggplot2)
library(dplyr)

# Load data from .csv files ----
fresh <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]

# Combine all data into a single data frame and change classifiers (species, organs, metabolites) into factors ----
allData <- rbind(fresh, frozenKR)
allData$species <- as.factor(allData$species)
allData$organ <- as.factor(allData$organ)
allData$metabolite <- as.factor(allData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame ----
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastafolia", "hastifolia"), collapse = '|')
excludeOrgans <- paste(c("flowers", "shoots"), collapse = '|')
# excludeMetabolites <- paste(c(), collapse = '|')
allData <- allData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) #>%>
# filter(!grepl(excludeMetabolites, metabolites))

# Fix naming errors ----
allData$species <- as.character(allData$species)
allData$species[allData$species=="RNA Seq"] <- "racemosa"
allData$species[allData$species=="havenesis"] <- "havanesis"
allData$species[allData$species=="hastafolia"] <- "hastifolia"
allData$species[allData$species=="havanesis"] <- "havanensis"
allData$species <- as.factor(allData$species)

# Average together any duplicate data points ----
allData <- allData %>%
  group_by(species, organ, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm))

# Define function to convert units of ppm to micromol/L ----
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

# Convert concentration and stError from ppm to mM for each data point ----
concentration_microM <- vector(mode="numeric", length=nrow(allData))
stError_microM <- vector(mode="numeric", length=nrow(allData))
for(i in 1:nrow(allData)){
  concentration_microM[i] <- ppm2microM(allData$concentration_ppm[i], allData$metabolite[i])
  stError_microM[i] <- ppm2microM(allData$stError_ppm[i], allData$metabolite[i])
}
allData$concentration_microM <- concentration_microM
allData$stError_microM <- stError_microM

# Set order of metabolites to appear in heatmaps based on pathway ----
allData$metabolite <- factor(allData$metabolite, levels=c(
  "chrysin", "chrysinG", "oroxylinA", "oroxyloside", "baicalein", "baicalin", "wogonin",
  "wogonoside", "acetoside", "apigenin", "apigeninG", "scutellarein", "scutellarin",
  "hispidulin", "hispidulinG")
) 

# Set colors to be used for metabolites across all plots ----
metaboliteColors <- c(
  "#8B0000", "#DC143C", "#FF7F50", "#FFD700", "#B8860B", "#BDB76B", "#808000", "#9ACD32", 
  "#2E8B57", "#66CDAA", "#2F4F4F", "#008080", "#4682B4", "#8A2BE2", "#8B008B")
names(metaboliteColors) <- levels(allData$metabolite)

# Define funciton to create metabolite bar plots ----
createPathwayPlots <- function(allData, metaboliteColors, speciesNames, metaboliteName){
  # Separate leaf and root data into separate objects, and make all root data negative
  leafData <- allData %>%
    filter(species %in% speciesNames & metabolite==metaboliteName & organ=="leaves")
  rootData <- allData %>%
    filter(species %in% speciesNames & metabolite==metaboliteName & organ=="roots") %>%
    mutate(concentration_microM=-concentration_microM) %>%
    mutate(stError_microM=-stError_microM)
  metaboliteData <- rbind(leafData, rootData)
  
  # Adjust order that species are displayed in plot to match order input with function call
  metaboliteData$species <- factor(metaboliteData$species, levels=speciesNames)
  
  # Calculate y-axis label values to use so that plot is symmetrical about x-axis
  roundUp <- function(x, accuracy){
    return(ceiling(x/ accuracy) * accuracy)
  }
  maxConc <- max(abs(metaboliteData$concentration_microM))
  if(maxConc<1){
    plotLimit <- roundUp(maxConc, 0.1)
  }else if(maxConc<10){
    plotLimit <- roundUp(maxConc, 1)
  }else if(maxConc<50){
    plotLimit <- roundUp(maxConc, 5)
  }else if(maxConc<100){
    plotLimit <- roundUp(maxConc, 10)
  }else{
    plotLimit <- roundUp(maxConc, 25)
  }
  
  breakValues <- c(-plotLimit, 0, plotLimit)
  
  print(maxConc)
  print(breakValues)
  
  barGraph <- ggplot(data=metaboliteData) +
    geom_col(mapping=aes(x=species, y=concentration_microM, fill=organ)) +
    geom_hline(yintercept=0, size=1) +
    geom_segment(x=0.4, y=breakValues[1]+breakValues[1]*0.005, xend=0.4, yend=breakValues[3]+breakValues[3]*0.005, size=1.5) +
    scale_fill_manual(values=c("#45abff", "#ff8745")) +
    scale_y_continuous(breaks = breakValues, labels = abs(breakValues), limits=c(breakValues[1], breakValues[3])) +
    theme_minimal() +
    theme(
      panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(),
      legend.position="none",
      axis.ticks.y=element_line(color="black", size=1),
      axis.ticks.length.y=unit(0.25, "cm"),
      axis.title=element_blank(),
      axis.text.x=element_blank(),
      axis.text=element_text(color="black", size=14)
    )
  
  return(barGraph)
}


testPlot <- createPathwayPlots(allData, metaboliteColors, c("baicalensis" ,"havanensis", "arenicola", "racemosa"), "oroxyloside")
print(testPlot)