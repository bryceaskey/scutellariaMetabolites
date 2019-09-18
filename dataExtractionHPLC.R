#install.packages("readxl")
library(readxl)

setwd("C:/Users/bca08_000/Documents/scutellaria/Test2 - Processing 2") # folder with HPLC excel files
filenames <- dir()

allData <- list()

for(i in 1:length(filenames)){
  fullSheet <- read_excel(filenames[i], sheet="Integration")
  startRow <- which(fullSheet[ , 1] == "Integration Results") + 4
  endRow <- nrow(fullSheet) - 1
  intResults <- fullSheet[startRow:endRow, ]
  colnames(intResults) <- paste(fullSheet[startRow - 3, ])
  allData[[i]] <- intResults
}

names(allData) <- filenames

exp1Data <- allData[-grep("-", names(allData))]
exp2Data <- allData[grep("-", names(allData))]

setwd("C:/Users/bca08_000/Documents/scutellaria") # folder containing excel file with dilution data
dilutionData <- read.csv("PeaksExp1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

exp1DataCorrected <- list()
exp1Count <- 1

for(i in 1:nrow(dilutionData)){
  sampleName <- paste(dilutionData$ï..Treatment[i], dilutionData$Replicate[i], sep="")
  index <- grep(sampleName, names(exp1Data))
  if(length(index) > 0){
    exp1DataCorrected[[i]] <- exp1Data[[index]]
    names(exp1DataCorrected)[i] <- names(exp1Data[index])
    if(dilutionData$Dilution[i] != ""){
      dilutionAmnt <- as.numeric(strsplit(dilutionData$Dilution[i], "%"))/100
      # only correcting area for dilution
      exp1DataCorrected[[i]]$Area <- as.numeric(exp1Data[[index]]$Area)/dilutionAmnt
    }
  }
}

# need to find sample with most rows - account for all possible metabolite peaks
peakCount <- vector(mode="numeric", length=length(allData))
for(i in 1:length(allData)){
  peakCount[i] <- nrow(allData[[i]])
}
allPeakNames <- allData[[which.max(peakCount)]]$'Peak Name'
allPeakNames[is.na(allPeakNames)] <- "Unknown"

# create new data frames to store height and area data
areaData <- data.frame(matrix(ncol=length(allPeakNames), nrow=length(exp1DataCorrected)))
heightData <- data.frame(matrix(ncol=length(allPeakNames), nrow=length(exp1DataCorrected)))

colnames(areaData) <- allPeakNames
colnames(heightData) <- allPeakNames

# pull out area data from all samples in experiment 1
for(i in 1:length(exp1DataCorrected)){
  sampleName <- names(exp1DataCorrected[i])
  peakNames <- exp1DataCorrected[[i]]$`Peak Name`
  peakAreas <- exp1DataCorrected[[i]]$Area
  peakHeights <- exp1DataCorrected[[i]]$Height 
  for(peak in peakNames){
    areaData[[peak]][i] <- peakAreas[which(peakNames == peak)]
    heightData[[peak]][i] <- peakHeights[which(peakNames == peak)]
  }
  rownames(areaData)[i] <- sampleName
  rownames(heightData)[i] <- sampleName
}

areaData[areaData=="n.a."] <- NA
heightData[areaData=="n.a."] <- NA

#example of plotting
barplot(as.numeric(areaData$IGS), names.arg=row.names(areaData), las=2)
barplot(as.numeric(heightData$IGS), names.arg=row.names(areaData), las=2)