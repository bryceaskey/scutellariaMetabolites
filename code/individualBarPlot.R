library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)

# Load data from .csv files
fresh <- read.csv("C:/Users/bca08_000/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/bca08_000/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]
wrightii <- read.csv("C:/Users/bca08_000/Documents/scutellariaMetabolites/data/preprocessed/20201007_wrightii.csv")[, 2:6]
suffrutescens <- read.csv("C:/Users/bca08_000/Documents/scutellariaMetabolites/data/preprocessed/20201119_suffrutescens.csv")[, 2:6]
cladeData <- read.csv("C:/Users/bca08_000/Documents/scutellariaMetabolites/data/phylo-tree-clades.csv")

# Remove barbata from fresh data - use only KR data
fresh <- fresh %>%
  filter(!grepl("barbata", species))

# Combine all data into a single data frame and change classifiers (species, organs, metabolites)
# into factors
allData <- rbind(fresh, frozenKR, wrightii)
allData$species <- as.factor(allData$species)
allData$organ <- as.factor(allData$organ)
allData$metabolite <- as.factor(allData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastafolia", "hastifolia", "arenicola", "havanesis"), collapse = '|')
excludeOrgans <- paste(c("flowers"), collapse = '|')
# excludeMetabolites <- paste(c(), collapse = '|')
allData <- allData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) #>%>
# filter(!grepl(excludeMetabolites, metabolites))

# Fix naming errors
allData$species <- as.character(allData$species)
allData$species[allData$species=="RNA Seq"] <- "racemosa"
allData$species[allData$species=="havenesis"] <- "havanensis"
allData$species[allData$species=="hastafolia"] <- "hastifolia"
allData$species[allData$species=="pekinesis"] <- "pekinensis var. alpina"
allData$species[allData$species=="indica"] <- "indica var. coccinea"
allData$species <- as.factor(allData$species)

allData$organ <- as.character(allData$organ)
allData$organ[allData$organ=="shoots"] <- "stems"
allData$organ <- as.factor(allData$organ)

# Drop unused factor levels
allData$species <- droplevels(allData$species)
allData$organ <- droplevels(allData$organ)

# Average together any duplicate data points
allData <- allData %>%
  group_by(species, organ, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm)) %>%
  ungroup()

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
  "#eff238", "#ffd320", "#ffb329", "#ff9040", "#ff6d5a", "#ff4b76", "#ff3291", "#e52dab", "#D12DE5")
names(metaboliteColors) <- levels(allData$metabolite)

createIndividualBars <- function(allData, metaboliteColors, indMetabolite){
  graphData <- allData %>%
    filter(metabolite==indMetabolite)
  
  graphData$species <- paste("S.", graphData$species)
  graphData$species <- factor(graphData$species, levels=c(
    "S. havanensis", "S. insignis", "S. indica var. coccinea", "S. barbata", "S. racemosa", "S. strigillosa", "S. dependens", "S. wrightii",
    #"S. arenicola", 
    "S. baicalensis", "S. tournefortii", "S. altissima", "S. leonardii", "S. pekinensis var. alpina"))
  
  indBarPlot <- ggplot(data=graphData, mapping=aes(x=species, y=concentration_microM, fill=organ)) +
    geom_col(position="dodge") +
    geom_errorbar(mapping=aes(ymin=concentration_microM-stError_microM, ymax=concentration_microM+stError_microM), color="black", width=0.2, position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#47acff", "#62c44d", "#ff8d4f"), name="Organ:", breaks=c("leaves", "stems", "roots"), labels=c("Leaf     ", "Stem     ", "Root")) +
    labs(y=paste(indMetabolite, "concentration (µmol/g FW)")) +
    theme_classic() +
    theme(panel.grid.major=element_line(size=0.5), panel.grid.minor=element_line(size=0.25),
          axis.title.x=element_blank(), axis.title.y=element_text(size=16), axis.text.y=element_text(size=14, color="black"),
          axis.text.x=element_text(size=14, face="italic", color="#000000", angle=90, hjust=1, vjust=0.5, margin=margin(25, 0, 0, 0)),
          legend.position="top", legend.title=element_text(size=14, face="bold"), legend.text=element_text(size=14))
  
  return(indBarPlot)
}

# Get clade data for all species to be plotted
allData$species <- factor(allData$species, levels=c(
  "havanensis", "insignis", "indica var. coccinea", "barbata", "racemosa", "strigillosa", "dependens", "wrightii",
  #"arenicola",
  "baicalensis", "tournefortii", "altissima", "leonardii", "pekinensis var. alpina")
)
speciesList <- vector(mode="character", length=length(levels(allData$species)))
cladeList <- vector(mode="numeric", length=length(levels(allData$species)))
for (i in 1:length(levels(allData$species))){
  species <- levels(allData$species)[i]
  speciesList[i] <- species
  clade <- cladeData$clade[cladeData$species==species]
  if (length(clade)>0 && !is.na(clade)){
    cladeList[i] <- clade
  }else{
    cladeList[i] <- NA
  }
} 
plotCladeData <- data.frame(x=1:length(levels(allData$species)), y=1, speciesList, cladeList)
plotCladeData$cladeList <- factor(plotCladeData$cladeList, levels=c(1, 2, 3, 4, 5))

# Create row of colored circles to represent phylogenetic clade
cladeLabels <- ggplot(data=plotCladeData) +
  geom_point(mapping=aes(x=x, y=y, fill=cladeList), shape=21, color="black", size=6) +
  scale_fill_manual(values=c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#9632B8"), drop=FALSE, na.value=NA) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=margin(-357,6,0,43.5,"pt"))


OroxylinA_plot <- createIndividualBars(allData, metaboliteColors, "OroxylinA")
OroxylinA_plot_withClades <- plot_grid(OroxylinA_plot, cladeLabels, nrow=2, ncol=1, rel_heights=c(1.5, 0.05))
ggsave(filename="C:/Users/bca08_000/Documents/scutellariaMetabolites/figures/indBarPlots/OroxylinA_plot.png",
       plot=OroxylinA_plot_withClades,
       device=png(),
       width=25, height=20, units="cm")

Oroxyloside_plot <- createIndividualBars(allData, metaboliteColors, "Oroxyloside")
Oroxyloside_plot_withClades <- plot_grid(Oroxyloside_plot, cladeLabels, nrow=2, ncol=1, rel_heights=c(1.5, 0.05))
ggsave(filename="C:/Users/bca08_000/Documents/scutellariaMetabolites/figures/indBarPlots/Oroxyloside_plot.png",
       plot=Oroxyloside_plot_withClades,
       device=png(),
       width=25, height=20, units="cm")