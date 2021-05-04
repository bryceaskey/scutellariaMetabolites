library(tidyverse)

# Load data from .csv files
fresh <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20190813_fresh.csv")[, 2:6]
frozenKR <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200117_frozenKR.csv")[, 2:6]
wrightii <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20201007_wrightii.csv")[, 2:6]
herbarium1_30 <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200214_herbarium1_30.csv")[, 2:6]
herbarium31_78 <- read.csv("C:/Users/Bryce/Documents/scutellariaMetabolites/data/preprocessed/20200812_herbarium31_78.csv")[, 2:6]

# Processing for fresh data ----
# Remove barbata from fresh data - use only KR data
fresh <- fresh %>%
  filter(!grepl("barbata", species))

# Combine all data into a single data frame and change classifiers (species, organs, metabolites)
# into factors
freshData <- rbind(fresh, frozenKR, wrightii)
freshData$species <- as.factor(freshData$species)
freshData$organ <- as.factor(freshData$organ)
freshData$metabolite <- as.factor(freshData$metabolite)

# Specify any species, organs, or metabolites to exclude, and remove from data frame
excludeSpecies <- paste(c("racemosa 071119", "racemosa MS", "racemosa SC", "hastafolia", "hastifolia"), collapse = '|')
excludeOrgans <- paste(c("flowers"), collapse = '|')
# excludeMetabolites <- paste(c(), collapse = '|')
freshData <- freshData %>%
  filter(!grepl(excludeSpecies, species)) %>%
  filter(!grepl(excludeOrgans, organ)) #>%>
# filter(!grepl(excludeMetabolites, metabolites))

# Fix naming errors
freshData$species <- as.character(freshData$species)
freshData$species[freshData$species=="RNA Seq"] <- "racemosa"
freshData$species[freshData$species=="havenesis"] <- "havanensis"
freshData$species[freshData$species=="hastafolia"] <- "hastifolia"
freshData$species[freshData$species=="pekinesis"] <- "pekinensis var. alpina"
freshData$species[freshData$species=="indica"] <- "indica var. coccinea"
freshData$species <- as.factor(freshData$species)

freshData$organ <- as.character(freshData$organ)
freshData$organ[freshData$organ=="shoots"] <- "stems"
freshData$organ <- as.factor(freshData$organ)

# Drop unused factor levels
freshData$species <- droplevels(freshData$species)
freshData$organ <- droplevels(freshData$organ)

# Average together any duplicate data points
freshData <- freshData %>%
  group_by(species, organ, metabolite) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm)) %>%
  ungroup()

# Adjust fresh ppm to correct for dilution
# Data is saved at 5000 ppm. Divide by 5 to calculate at 1000 ppm (= umol/1 g FW)
freshData <- freshData %>%
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
concentration_microM <- vector(mode="numeric", length=nrow(freshData))
stError_microM <- vector(mode="numeric", length=nrow(freshData))
for(i in 1:nrow(freshData)){
  concentration_microM[i] <- ppm2microM(freshData$concentration_ppm[i], freshData$metabolite[i])
  stError_microM[i] <- ppm2microM(freshData$stError_ppm[i], freshData$metabolite[i])
}
freshData$concentration_microM <- concentration_microM
freshData$stError_microM <- stError_microM

# Filter data to only include necessary data
freshData <- freshData[, c(1:3, 6:7)]
freshData <- freshData %>%
  pivot_wider(names_from=metabolite, values_from=c(concentration_microM, stError_microM))
#write.csv(freshData, file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/metaboliteConcentrations.csv")

# Processing for herbarium samples ----
herbariumData <- rbind(herbarium1_30, herbarium31_78)
herbariumData$species <- factor(herbariumData$species)
herbariumData$organ <- factor(herbariumData$organ)
herbariumData$metabolite <- factor(herbariumData$metabolite)

excludeOrgans <- c("roots")
herbariumData <- herbariumData %>%
  filter(!organ %in% excludeOrgans)
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

herbariumData$organ <- as.character(herbariumData$organ)
herbariumData$organ[herbariumData$organ=="shoots"] <- "stems"
herbariumData$organ <- as.factor(herbariumData$organ)

# Average together duplicate species, and leaf and shoot data
herbariumData <- herbariumData %>%
  group_by(species, metabolite, organ) %>%
  summarise(concentration_ppm=mean(concentration_ppm), stError_ppm=mean(stError_ppm)) %>%
  ungroup()

# Adjust herbarium ppm to correct for dilution
# Data is saved at 1000 ppm. Divide by 10 to calculate at 100 ppm (= umol/0.1 g DW)
herbariumData <- herbariumData %>%
  transmute(
    species=species,
    organ=organ,
    metabolite=metabolite,
    concentration_ppm=concentration_ppm/10,
    stError_ppm=stError_ppm/10
  )

# Convert concentration and stError from ppm to mM for each data point
concentration_microM <- vector(mode="numeric", length=nrow(herbariumData))
stError_microM <- vector(mode="numeric", length=nrow(herbariumData))
for(i in 1:nrow(herbariumData)){
  concentration_microM[i] <- ppm2microM(herbariumData$concentration_ppm[i], herbariumData$metabolite[i])
  stError_microM[i] <- ppm2microM(herbariumData$stError_ppm[i], herbariumData$metabolite[i])
}
herbariumData$concentration_microM <- concentration_microM
herbariumData$stError_microM <- stError_microM

# Filter data to only include necessary data
herbariumData <- herbariumData[, c(1:3, 6)]
herbariumData <- herbariumData %>%
  pivot_wider(names_from=metabolite, values_from=c(concentration_microM))
write.csv(herbariumData, file="C:/Users/Bryce/Documents/scutellariaMetabolites/data/metaboliteConcentrations.csv")