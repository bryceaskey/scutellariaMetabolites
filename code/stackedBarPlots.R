# TODO: if a species doesn't produce any of the 15 metabolites in an organ, the bars will be
# shifted left but labels won't be. Add placeholder data point and then photoshop out later?

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
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastafolia", "hastifolia"), collapse = '|')
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
allData$species[allData$species=="pekinesis"] <- "pekinensis"
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

# Set order of metabolites to appear in heatmaps based on pathway
allData$metabolite <- factor(allData$metabolite, levels=c(
  "chrysin", "chrysinG", "baicalein", "baicalin", "oroxylinA", "oroxyloside", "wogonin",
  "wogonoside", "acetoside", "apigenin", "apigeninG", "scutellarein", "scutellarin",
  "hispidulin", "hispidulinG")
) 
# Set order of metabolites to appear in heatmaps based on heirarchical clustering
#allData$metabolite <- factor(allData$metabolite, levels=c(
#  "acetoside", "hispidulinG", "baicalin", "chrysinG", "oroxylinA", "oroxyloside", "hispidulin",
#  "baicalein", "wogonin", "wogonoside", "scutellarein", "scutellarin", "chrysin", "apigenin",
#  "apigeninG"
#))

# Create new column w/ factor nums - for bar plot section labeling
allData$metNum <- as.numeric(allData$metabolite) 

# Set colors to be used for metabolites across all plots
metaboliteColors <- c(
  "#e50000", "#e65235", "#da6900", "#f2a155", "#CF9C00", "#C9C300", "#83c400", "#2dbf00", 
  "#00D8E5", "#00b3e5", "#4232f0", "#1000D2", "#a032f0", "#a200bf", "#bf009c")
names(metaboliteColors) <- levels(allData$metabolite)

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
  speciesList <- levels(allData$species)
  
  organData <- allData[order(allData$metNum), ]
  organData <- organData %>%
    filter(grepl(plantOrgan, organ)) %>%
    filter(concentration_microM > 0) %>%
    group_by(species) %>%
    mutate(text_y = sum(concentration_microM) - (cumsum(concentration_microM) - concentration_microM/2))
  organData$species <- paste("S.", organData$species)
  organData$species <- factor(organData$species, levels=c(
    "S. baicalensis", "S. strigillosa", "S. dependens", "S. indica", "S. barbata", "S. insignis", "S. racemosa", 
    "S. arenicola", "S. havanensis", "S. altissima", "S. tournefortii", "S. leonardii", "S. pekinensis")
  )
  
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
    ylim(0, 400) +
    labs(title=paste("Flavonoid concentrations in Scutellaria", plantOrgan),
         x="Species", 
         y=expression(paste("Concentration (", mu, "M)", sep=""))) +
    scale_fill_manual(values=metaboliteColors) +
    geom_text_repel(mapping=aes(label=metNum, x=text_x, y=text_y), hjust=1, direction="y", nudge_x=-0.2) +
    if(plantOrgan=="roots"){
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_text(size=20, color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(30, 0, 0, 0)),
            axis.text.y=element_text(color="#000000"),
            panel.background=element_rect(fill="#ffe0cf"),
            text=element_text(size=20))
    }else if(plantOrgan=="shoots"){
      theme(legend.position="none", 
            axis.title.x=element_blank(),
            axis.text.x=element_text(size=20, color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(30, 0, 0, 0)),
            axis.text.y=element_text(color="#000000"),
            panel.background=element_rect(fill="#d5ffcc"),
            text=element_text(size=20))
    }else{
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_text(size=20, color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(30, 0, 0, 0)),
            axis.text.y=element_text(color="#000000"),
            panel.background = element_rect(fill="#cce8ff"),
            text=element_text(size=20))
    }
  print(chart)
}

rootPlot <- createStackedBars(allData, metaboliteColors, "roots")
shootPlot <- createStackedBars(allData, metaboliteColors, "shoots")
leafPlot <- createStackedBars(allData, metaboliteColors, "leaves")

#justData <- plot_grid(leafPlot, shootPlot, rootPlot, nrow=3, ncol=1, rel_heights = c(1, 1, 1.1))
legend <- createLegend(allData, metaboliteColors, legendOrientation="vertical")

finalRootFigure <- plot_grid(rootPlot, legend,
  nrow=1, ncol=2, 
  rel_widths=c(1, 0.15))
print(finalRootFigure)

finalShootFigure <- plot_grid(shootPlot, legend,
  nrow=1, ncol=2, 
  rel_widths=c(1, 0.15))
print(finalShootFigure)

finalLeafFigure <- plot_grid(leafPlot, legend,
  nrow=1, ncol=2,
  rel_widths=c(1, 0.15))
print(finalLeafFigure)

# Export at size 1800x1000