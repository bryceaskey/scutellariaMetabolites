# Function to identify metabolites in HPLC chromatogram data.
# Associates peak retention times from standards with peak retention times from mixture.
# Generate accuracy statistics (r^2, RMSE) to allow for identification of most correct association.

# load data from .csv file
standards <- read.csv("metaboliteRetentionTimes.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, skip=1)[1:15, 1:2]
unknowns <- read.csv("metaboliteRetentionTimes.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, skip=1)[1:15, 5:6]

for(i in 1:1000){
  unknownRTs <- unknowns[sample(nrow(unknowns)), ] # shuffle row-wise
  standardRTs <- standards[, 2]
  matchingRTs <- vector(length=15)
  matchingMix <- vector(length=15)
  for(j in 1:15){
    unknownRT = as.numeric(unknowns[j, 2])
    metaboliteMatch <- which.min(abs(as.numeric(standardRTs) - unknownRT))
    standardRTs[metaboliteMatch] <- 1000 # prevent reconsideration
    matchingRTs[j] <- unknownRTs[metaboliteMatch, 2]
    matchingMix[j] <- unknownRTs[metaboliteMatch, 1]
  }
  outputMatches <- data.frame(standards[, 1:2], matchingRTs, matchingMix)
}