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
  # Create color legend for metabolite identification
  x <- ggplot(filter(allData, species==levels(allData$species)[1])) +
    geom_bar(mapping=aes(x="", y=meanConc, fill=metabolite), stat="identity") +
    theme(legend.position="bottom", legend.direction="horizontal",
          legend.text=element_text(size=16), legend.title=element_text(size=20)) +
    labs(fill="Flavonoid") +
    scale_fill_manual(values=metaboliteColors, 
                      labels=paste(seq(1, 15), ". ", sapply(names(metaboliteColors), capStr), sep=""))
  return(get_legend(x))
}

# Function to create stacked bar plot from metabolite data ----
createPlot <- function(filteredData, legend){
  # Calculate x and y locations for metabolite labels
  filteredData <- filteredData %>%
    group_by(species) %>%
    arrange(metNum, .by_group=TRUE) %>%
    mutate(text_y = sum(meanConc) - (cumsum(meanConc) - meanConc/2)) %>%
    mutate(text_x = as.numeric(species) - 0.25)
  filteredData$metNum[filteredData$meanConc == 0] = ""
  
  # Create stacked bar plot
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
  finalFigure <- plot_grid(chart, legend, nrow=2, ncol=1,rel_heights=c(1, 0.20))
  return(finalFigure)
  #TODO: add in legend
}
