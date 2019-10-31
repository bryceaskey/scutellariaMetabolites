# install any necessary packages if they aren't already installed
packageList <- c("ggplot2", "dplyr", "tidyr", "magrittr", "gridExtra", "readxl")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0){
  install.packages(newPackages)
}

# load packages into workspace
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(gridExtra)
library(readxl)

# set working directory to folder containing excel file
setwd("C:/Users/bca08_000/Documents/scutellaria")

# save excel file as a dataframe
fullSheet <- read_excel("HPLC Ru Data-2.xlsx", .name_repair="minimal")

# extract peak sums from dataframe into individual vectors
allSums <- as.numeric(fullSheet[["SUM"]])
wtSums <- allSums[2:7]
B2Ox1Sums <- allSums[8:13]
B2Ox2Sums <- allSums[14:19]
B2Ox3Sums <- allSums[20:25]


print(t.test(wtSums, B2Ox1Sums))
print(t.test(wtSums, B2Ox2Sums))
print(t.test(wtSums, B2Ox3Sums))