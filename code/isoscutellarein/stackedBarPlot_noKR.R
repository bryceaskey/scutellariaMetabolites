library(tidyverse)
library(ggrepel)
library(cowplot)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20190813_fresh.csv")
wrightii <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20210119_wrightii.csv")
suffrutescens <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201119_suffrutescens.csv")

# Combine all data into a single data frame
allData <- rbind(fresh, wrightii, suffrutescens)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa_071119", "racemosa_MS", "racemosa_SC", "hastifolia", "havanensis", "arenicola", "suffrutescens"), collapse="|")
excludeOrgans <- paste(c("flowers"), collapse="|")
excludeMetabolites <- paste(c("acteoside", "isoscutellarin"), collapse="|")
allData <- allData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) %>%
  filter(!grepl(excludeMetabolites, metabolite))

# Fix naming errors
allData$species[allData$species=="racemosa_RNAseq"] <- "racemosa"

# Change classifiers (species, organs, metabolites) into factors
allData$species <- factor(allData$species)
allData$organ <- factor(allData$organ)
allData$metabolite <- factor(allData$metabolite)

# Adjust fresh ppm to correct for dilution
# Data is saved at 5000 ppm. Divide by 5 to calculate at 1000 ppm (= umol/1 g FW)
allData <- allData %>%
  transmute(
    species=species,
    organ=organ,
    metabolite=metabolite,
    concentration_ppm=concentration_ppm/5,
    stError_ppm=stError_ppm/5
  )

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
    }else if(metaboliteName=="isoscutellarin"){
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

allData$metabolite <- as.character(allData$metabolite)
allData$metabolite[allData$metabolite=="apigeninG"] <- "apigenin 7-G"
allData$metabolite[allData$metabolite=="chrysinG"] <- "chrysin 7-G"
allData$metabolite[allData$metabolite=="hispidulinG"] <- "hispiduloside"
allData$metabolite[allData$metabolite=="oroxylinA"] <- "oroxylin A"
allData$metabolite <- factor(allData$metabolite)

# Capitalize first letter of each flavonoid name
capString <- function(string) {
  c <- strsplit(string, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2), sep="", collapse=" ")
}
allData$metabolite <- as.character(allData$metabolite)
allData$metabolite <- sapply(allData$metabolite, capString)

# Set order of metabolites to appear in heatmaps based on pathway
allData$metabolite <- factor(allData$metabolite, levels=c(
  "Apigenin", "Apigenin 7-G", "Scutellarein", "Scutellarin", "Hispidulin", "Hispiduloside",
  "Chrysin", "Chrysin 7-G", "Baicalein", "Baicalin", "Oroxylin A", "Oroxyloside", "Wogonin", "Wogonoside")
) 

# Create new column w/ factor nums - for bar plot section labeling
allData$metNum <- as.numeric(allData$metabolite) 

# Set colors to be used for metabolites across all plots
metaboliteColors <- c(
  "#4726dd", "#006bff", "#008afe", "#009ec2", "#00ad76", "#169E18", 
  "#eff238", "#ffd320", "#ffb329", "#ff9040", "#ff6d5a", "#ff4b76", "#ff3291", "#e52dab")
names(metaboliteColors) <- levels(allData$metabolite)

# Function to create color legend for scaled pie charts
createLegend <- function(allData, metaboliteColors, legendOrientation="horizontal"){
  x <- ggplot(filter(allData, species==levels(allData$species)[1])) +
    geom_bar(mapping=aes(x="", y=concentration_microM, fill=metabolite), stat="identity") +
    theme(legend.position="bottom", 
          legend.direction=legendOrientation, 
          legend.text=element_text(size=8), 
          legend.title=element_text(size=8, face="bold"),
          legend.key.size=unit(0.5, "cm"),
          legend.box.margin=unit(c(-2.5,0,0,0.1), "cm")) +
    labs(fill="Flavonoid") +
    scale_fill_manual(name="Flavone:", values=metaboliteColors, labels=paste(seq(1, 14), ". ", names(metaboliteColors), sep=""))
  legend <- get_legend(x)
  return(legend)
}

# Function to create stacked bar charts
createStackedBars <- function(allData, metaboliteColors, plantOrgan){
  organData <- allData[order(allData$metNum), ]
  organData <- organData %>%
    filter(grepl(plantOrgan, organ)) %>%
    filter(concentration_microM > 0) %>%
    group_by(species) %>%
    mutate(text_y = sum(concentration_microM) - (cumsum(concentration_microM) - concentration_microM/2))
  
  if(length(levels(droplevels(organData$species))) != length(levels(organData$species))){
    print("TRUE")
    for(speciesName in levels(organData$species)){
      if(sum(grepl(speciesName, levels(droplevels(organData$species)))) == 0){
        print("TRUE")
        NA_df <- data.frame(species=speciesName, organ=NA, metabolite=NA, concentration_ppm=0,
                            stError_ppm=NA, concentration_microM=0, stError_microM=NA,
                            metNum=NA, text_x=NA, text_y=NA)
        organData <- rbind(organData, NA_df)
        organData$species <- factor(organData$species)
      }
    }
  }
  
  organData$species <- paste("S.", organData$species)
  organData$species <- str_wrap(organData$species, width=16)
  organData$species <- factor(organData$species)
  
  organData <- organData %>%
    group_by(species) %>%
    mutate(text_x = as.numeric(species) - 0.25)
  
  levels(organData$species) <- str_wrap(levels(organData$species), width=16)
  
  chart <- ggplot(data=organData, mapping=aes(x=species, y=concentration_microM, fill=metabolite)) +
    geom_bar(position="stack", stat="identity", width=0.45) +
    #ylim(0, 400) +
    labs(x="Species",
         y=paste("Concentration in ", plantOrgan, " (µmol/g FW)",  sep="")) +
    scale_fill_manual(values=metaboliteColors) +
    theme(axis.text.x=element_text(face="italic")) +
    geom_text_repel(mapping=aes(label=metNum, x=text_x, y=text_y), nudge_x=-0.25, direction="y", size=6/.pt,
                    segment.size=0.2, box.padding=0.15, segment.alpha=0.6, min.segment.length=0.25) +
    if(plantOrgan=="roots"){
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_text(size=8, color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(4, 0, 0, 0)),
            axis.text.y=element_text(color="#000000", size=8), axis.title.y=element_text(color="#000000", size=8),
            panel.background=element_rect(fill="#ffe0cf"),
            plot.margin=margin(0, 0, 0, 0, unit="cm"))
    }else if(plantOrgan=="stems"){
      theme(legend.position="none", 
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            #axis.text.x=element_text(size=18, color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(25, 0, 0, 0)),
            axis.text.y=element_text(color="#000000", size=8), axis.title.y=element_text(color="#000000", size=8),
            panel.background=element_rect(fill="#d5ffcc"),
            plot.margin=margin(0, 0, 0.15, 0, unit="cm"))
    }else{
      theme(legend.position="none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            #axis.text.x=element_text(size=18, color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(25, 0, 0, 0)),
            axis.text.y=element_text(color="#000000", size=8), axis.title.y=element_text(color="#000000", size=8),
            panel.background = element_rect(fill="#cce8ff"),
            plot.margin=margin(0, 0, 0.15, 0, unit="cm"))
    }
}

rootPlot <- createStackedBars(allData, metaboliteColors, "roots")
shootPlot <- createStackedBars(allData, metaboliteColors, "stems")
leafPlot <- createStackedBars(allData, metaboliteColors, "leaves")


legend <- createLegend(allData, metaboliteColors, legendOrientation="vertical")

allOrganPlot <- plot_grid(leafPlot, shootPlot, rootPlot, nrow=3, ncol=1, rel_heights=c(1,1,1.3), align="v", axis="l")
completePlot <- plot_grid(allOrganPlot, legend, nrow=1, ncol=2, rel_widths=c(1,0.3))

ggsave(filename="C:/Users/Bryce/Research/scutellariaMetabolites/figures/0-isoscutellarin/Figure_2.pdf",
      plot=completePlot,
      device=pdf(),
      width=5.1, height=9, units="in")

# Generate reader-friendly table of data
tableConcData <- allData %>% 
  pivot_wider(id_cols=c("species", "organ"), names_from=c("metabolite"), values_from=c("concentration_microM"))

tableErrorData <- allData %>%
  pivot_wider(id_cols=c("species", "organ"), names_from=c("metabolite"), values_from=c("stError_microM"))

tableData <- tableConcData
tableData[ , 3:ncol(tableConcData)] <- NA
tableData$species <- paste("S.", tableData$species)
tableData$organ <- as.character(tableData$organ)
tableData$organ <- sapply(tableData$organ, capString)

for(metabolite in colnames(tableConcData)[3:ncol(tableConcData)]){
  concData <- format(round(tableConcData[[metabolite]], 2), nsmall = 2)
  errorData <- format(round(tableErrorData[[metabolite]], 2), nsmall = 2)
  tableData[[metabolite]] <- paste(concData, "\u00B1", errorData)
}

write.csv(tableData, file="C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/tables/noKR.csv", row.names=FALSE)
