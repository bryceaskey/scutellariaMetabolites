library(tidyverse)
library(ggrepel)
library(cowplot)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20190813_fresh.csv")
wrightii <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201007_wrightii.csv")
suffrutescens <- read.csv("C:/Users/Bryce/Research/scutellariaMetabolites/data/hplc/preprocessed/20201119_suffrutescens.csv")

# Combine all data into a single data frame
allData <- rbind(fresh, wrightii, suffrutescens)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa_071119", "racemosa_MS", "racemosa_SC", "hastifolia", "havanensis"), collapse="|")
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

# Set order of organs to appear in plot
allData$organ <- factor(allData$organ, levels=c("leaves", "stems", "roots"))

# Set colors to be used for metabolites across all plots
metaboliteColors <- c(
  "#4726dd", "#006bff", "#008afe", "#009ec2", "#00ad76", "#169E18", 
  "#eff238", "#ffd320", "#ffb329", "#ff9040", "#ff6d5a", "#ff4b76", "#ff3291", "#e52dab")
names(metaboliteColors) <- levels(allData$metabolite)

createIndividualBars <- function(allData, metaboliteColors, indMetabolite, axisLabels=TRUE, legend=TRUE){
  graphData <- allData %>%
    filter(metabolite==indMetabolite)
  
  graphData$species <- paste("S.", graphData$species)
  graphData$species <- str_wrap(graphData$species, width=16)
  graphData$species <- factor(graphData$species)
  
  indBarPlot <- ggplot(data=graphData, mapping=aes(x=species, y=concentration_microM, fill=organ)) +
    geom_col(position="dodge") +
    geom_errorbar(mapping=aes(ymin=concentration_microM-stError_microM, ymax=concentration_microM+stError_microM), color="black", width=0.35, size=0.35, position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#47acff", "#62c44d", "#ff8d4f"), name="Organ:", breaks=c("leaves", "stems", "roots"), labels=c("Leaf     ", "Stem     ", "Root")) +
    labs(y=paste(indMetabolite, "concentration (µmol/g FW)")) +
    theme_classic() +
    theme(panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
          axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text.y=element_text(size=8, color="black"),
          axis.text.x=element_text(size=8, face="italic", color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(4, 0, 0, 0)),
          legend.position="top", legend.title=element_text(size=8), legend.text=element_text(size=8), legend.key.size=unit(0.35, "cm"), legend.box.margin=margin(-5,0,-8,0),
          axis.ticks=element_line(color="black")
    )
  
  if(axisLabels == FALSE){
    indBarPlot <- indBarPlot + theme(axis.text.x=element_blank())
  }
  
  if(legend == FALSE){
    indBarPlot <- indBarPlot + theme(legend.position="none")
  }
  
  
  return(indBarPlot)
}

OroxylinAPlot <- createIndividualBars(allData, metaboliteColors, "Oroxylin A", axisLabels=FALSE, legend=TRUE)

OroxylosidePlot <- createIndividualBars(allData, metaboliteColors, "Oroxyloside", axisLabels=TRUE, legend=FALSE)

combinedPlot <- plot_grid(OroxylinAPlot, OroxylosidePlot, nrow=2, ncol=1, rel_heights=c(1,1.35))
ggsave(filename="C:/Users/Bryce/Research/scutellariaMetabolites/figures/indBarPlots/combinedPlot_noKR.pdf",
       plot=combinedPlot,
       device=pdf(),
       width=5, height=6, units="in")