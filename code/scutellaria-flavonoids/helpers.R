# Function to capitalize the first letter of an input string
# Used by createPlot function to capitalize first letter of metabolite names
capStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

# List metabolites to match order of metNum in allData
metaboliteOrder <- c("oroxyloside", "oroxylinA", "hispidulinG", "hispidulin", "chrysin",
                     "chrysinG", "apigenin", "apigeninG", "acetoside", "scutellarein",
                     "scutellarin", "baicalin", "baicalein", "wogonin", "wogonoside")

# Set colors to be used for metabolites across all plots
metaboliteColors <- c("#8B0000", "#DC143C", "#FF7F50", "#FFD700", "#B8860B", "#BDB76B",
                      "#808000", "#9ACD32", "#2E8B57", "#66CDAA", "#2F4F4F", "#008080",
                      "#4682B4", "#8A2BE2", "#8B008B")

# Create dataframe of metabolite names and colors to be used in legend creation
names(metaboliteColors) <- metaboliteOrder

# Function to create color legend for metabolite identification
createLegend <- function(allData, metaboliteColors){
  x <- ggplot(filter(allData, species==levels(allData$species)[1])) +
    geom_bar(mapping=aes(x="", y=meanConc, fill=metabolite), stat="identity") +
    theme(legend.position="bottom", legend.direction="vertical",
          legend.text=element_text(size=16), legend.title=element_text(size=20)) +
    labs(fill="Flavonoid") +
    scale_fill_manual(values=metaboliteColors, 
                      labels=paste(seq(1, 15), ". ", names(metaboliteColors), sep=""))
  legend <- plot_grid(get_legend(x), nrow=1, ncol=1)
  return(legend)
}

# Function to filter data based on ui inputs, and calculate data and label positions ----
filterData <- function(allData, selectedOrgan, selectedMetabolites){
  filteredData <- allData %>%
    filter(grepl(selectedOrgan, organ)) %>%
    filter(grepl(paste(paste("^", selectedMetabolites, "$", sep=""), collapse="|"), metabolite)) %>%
    group_by(species) %>%
      arrange(metNum, .by_group=TRUE) %>%
      mutate(text_y = sum(meanConc) - (cumsum(meanConc) - meanConc/2)) %>%
      mutate(text_x = as.numeric(species) - 0.25)
  filteredData$metNum[filteredData$meanConc == 0] <- ""
  return(filteredData)
}

# Function to create stacked bar plot from metabolite data ----
createPlot <- function(filteredData, metaboliteColors){
  chart <- ggplot(data=filteredData, mapping=aes(x=species, y=meanConc, fill=metabolite)) +
    geom_bar(position="stack", stat="identity", width=0.5) +
    labs(x="Species", y="Concentration (ppm)") +
    scale_fill_manual(values=metaboliteColors) +
    geom_text_repel(mapping=aes(label=metNum, x=text_x, y=text_y), hjust=1, direction="y", nudge_x=-0.2) +
    if(filteredData$organ[[1]]=="Roots"){
      theme(legend.position="none",
            panel.background = element_rect(fill="#ffe0cf"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=20))
    }else if(filteredData$organ[[1]]=="Shoots"){
      theme(legend.position="none", 
            panel.background = element_rect(fill="#d5ffcc"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=20))
    }else{
      theme(legend.position="none",
            panel.background = element_rect(fill="#cce8ff"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=20))
    }
  return(chart)
}

# Function to create an empty plot to display before user has selected any flavonoids ----
createEmptyPlot <- function(allData, selectedOrgan){
  chart <- ggplot(data=allData, mapping=aes(x=species, y=0)) +
    labs(x="Species", y="Concentration (ppm)") + 
    ylim(0, 1) +
    if(selectedOrgan == "Roots"){
      theme(legend.position="none",
            panel.background = element_rect(fill="#ffe0cf"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=20))
    }else if(selectedOrgan == "Shoots"){
      theme(legend.position="none",
            panel.background = element_rect(fill="#d5ffcc"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=20))
    }else{
      theme(legend.position="none",
            panel.background = element_rect(fill="#cce8ff"),
            panel.grid.minor.y=element_blank(), panel.grid.major.x=element_blank(),
            text=element_text(size=20))
    }
  return(chart)
}

# Function to refine data for display
refineData <- function(allData, selectedOrgan, selectedMetabolites, selectedSpecies){
  filteredData <- allData %>%
    filter(grepl(selectedOrgan, organ)) %>%
    filter(grepl(paste(paste("^", selectedMetabolites, "$", sep=""), collapse="|"), metabolite)) %>%
    group_by(species) %>%
      arrange(metNum, .by_group=TRUE) %>%
      mutate(text_y = sum(meanConc) - (cumsum(meanConc) - meanConc/2)) %>%
      mutate(text_x = as.numeric(species) - 0.25) %>%
    filter(grepl(selectedSpecies, species))
  refinedData <- filteredData[, c(8, 7, 3, 4)]
  refinedData$metNum <- as.character(refinedData$metNum)
  colnames(refinedData) <- c("Flavonoid #", "Flavonoid name", "Concentration (ppm)", "Standard error")
  return(refinedData)
}
