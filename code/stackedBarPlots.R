library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]

# Combine all data into a single data frame and change classifiers (species, organs, metabolites)
# into factors
allData <- rbind(fresh, frozenKR)
allData$species <- as.factor(allData$species)
allData$organ <- as.factor(allData$organ)
allData$metabolite <- as.factor(allData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC"), collapse = '|')
excludeOrgans <- paste(c("flowers"), collapse = '|')
# excludeMetabolites <- paste(c(), collapse = '|')
allData <- allData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) #>%>
# filter(!grepl(excludeMetabolites, metabolites))

# Fix naming errors
allData$species <- as.character(allData$species)
allData$species[allData$species=="RNA Seq"] <- "racemosa"
allData$species[allData$species=="havenesis"] <- "havanesis"
allData$species[allData$species=="hastafolia"] <- "hastifolia"
allData$species[allData$species=="havanesis"] <- "havanensis"
allData$species <- as.factor(allData$species)

# Average together any duplicate data points
allData <- allData %>%
  group_by(species, organ, metabolite) %>%
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
concentration_microM <- vector(mode="numeric", length=nrow(allData))
stError_microM <- vector(mode="numeric", length=nrow(allData))
for(i in 1:nrow(allData)){
  concentration_microM[i] <- ppm2microM(allData$concentration_ppm[i], allData$metabolite[i])
  stError_microM[i] <- ppm2microM(allData$stError_ppm[i], allData$metabolite[i])
}
allData$concentration_microM <- concentration_microM
allData$stError_microM <- stError_microM

# Set order of species and metabolites to appear in heatmaps
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin",
                     "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein",
                     "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")
levels(allData$metabolite) <- metaboliteOrder

# Create new column w/ factor nums - for bar plot section labeling
allData$metNum <- as.numeric(allData$metabolite) 

# Set colors to be used for metabolites across all plots
metaboliteColors <- c("#8B0000", "#DC143C", "#FF7F50", "#FFD700", "#B8860B", "#BDB76B", "#808000", "#9ACD32", "#2E8B57", "#66CDAA", "#2F4F4F", "#008080", "#4682B4", "#8A2BE2", "#8B008B")
names(metaboliteColors) <- metaboliteOrder

# Function to create color legend for scaled pie charts
createLegend <- function(allData, metaboliteColors, legendOrientation="horizontal"){
  x <- ggplot(filter(allData, species==levels(allData$species)[1])) +
    geom_bar(mapping=aes(x="", y=concentration_microM, fill=metabolite), stat="identity") +
    theme(legend.position="bottom", legend.direction=legendOrientation, legend.text=element_text(size=12), legend.title=element_text(size=16)) +
    labs(fill="Flavonoid") +
    scale_fill_manual(values=metaboliteColors, labels=paste(seq(1, 15), ". ", names(metaboliteColors), sep=""))
  legend <- get_legend(x)
  return(legend)
}

# Function to create stacked bar charts
createStackedBars <- function(allData, metaboliteColors, plantOrgan){
  allData$species <- as.character(allData$species)
  allData$species <- factor(allData$species, levels=c("baicalensis", "havanensis", "arenicola", "hastifolia", 
  "dependens", "strigillosa", "barbata", "indica", "insignis", "racemosa", "tournefortii", "altissima", "leonardii", "pekinesis"))
  speciesList <- levels(allData$species)
  
  print(levels(allData$species))
  
  organData <- allData[order(allData$metNum), ]
  organData <- organData %>%
    filter(grepl(plantOrgan, organ)) %>%
    filter(concentration_microM > 0) %>%
    group_by(species) %>%
    mutate(text_y = sum(concentration_microM) - (cumsum(concentration_microM) - concentration_microM/2))
  organData$species <- factor(organData$species)
  
  if(length(speciesList) != length(levels(organData$species))){
    for(speciesName in speciesList){
      if(sum(grepl(speciesName, levels(organData$species))) == 0){
        NA_df <- data.frame(species=speciesName, organ=NA, metabolite=NA, concentration_ppm=0,
                            stError_ppm=NA, concentration_microM=0, stError_microM=NA,
                            metNum=NA, text_x=NA, text_y=NA)
        organData <- rbind(organData, NA_df)
        organData$species <- as.factor(organData$species)
      }
    }
  }
  
  
  organData <- organData %>%
    group_by(species) %>%
    mutate(text_x = as.numeric(species) - 0.25)
  print(organData)
  
  chart <- ggplot(data=organData, mapping=aes(x=species, y=concentration_microM, fill=metabolite)) +
    geom_bar(position="stack", stat="identity", width=0.5) +
    # ylim(0, 100) +
    labs(title=paste("Flavonoid concentrations in Scutellaria", plantOrgan),
         x="Species", 
         y=expression(paste("Concentration (", mu, "M)", sep=""))) +
    scale_fill_manual(values=metaboliteColors) +
    geom_text_repel(mapping=aes(label=metNum, x=text_x, y=text_y), hjust=1, direction="y", nudge_x=-0.2) +
    if(plantOrgan=="roots"){
      theme(legend.position="none",
            axis.text.x=element_text(color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(30, 0, 0, 0)),
            axis.text.y=element_text(color="#000000"),
            axis.title.x=element_blank(),
            panel.background=element_rect(fill="#ffe0cf"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=18))
    }else if(plantOrgan=="shoots"){
      theme(legend.position="none", 
            #axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_text(color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(30, 0, 0, 0)),
            axis.text.y=element_text(color="#000000"),
            axis.title.x=element_blank(),
            panel.background=element_rect(fill="#d5ffcc"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=18))
    }else{
      theme(legend.position="none",
            #axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_text(color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(30, 0, 0, 0)),
            axis.text.y=element_text(color="#000000"),
            axis.title.x=element_blank(),
            panel.background = element_rect(fill="#cce8ff"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=18))
    }
}

rootPlot <- createStackedBars(allData, metaboliteColors, "roots")
shootPlot <- createStackedBars(allData, metaboliteColors, "shoots")
leafPlot <- createStackedBars(allData, metaboliteColors, "leaves")

#justData <- plot_grid(leafPlot, shootPlot, rootPlot, nrow=3, ncol=1, rel_heights = c(1, 1, 1.1))
legend <- createLegend(allData, metaboliteColors, legendOrientation="vertical")

finalRootFigure <- plot_grid(rootPlot, legend,
  nrow=1, ncol=2, 
  rel_widths=c(1, 0.20))
print(finalRootFigure)

finalShootFigure <- plot_grid(shootPlot, legend,
  nrow=1, ncol=2, 
  rel_widths=c(1, 0.20))
print(finalShootFigure)

finalLeafFigure <- plot_grid(leafPlot, legend,
  nrow=1, ncol=2,
  rel_widths=c(1, 0.20))
print(finalLeafFigure)