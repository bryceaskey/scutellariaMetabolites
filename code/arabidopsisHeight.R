library(ggplot2)
library(plyr)

heightData <- read.csv(file="C:/Users/bca08_000/Documents/scutellaria/data/Col and A2 -1 combined.csv")
colnames(heightData)  <- c("genotype", "group", "height")
heightData[, 2] <- as.factor(heightData[, 2])

# To separate height data by genotype only --------------------------------------------------------
meanHeights <- ddply(heightData, c("genotype"), summarise, mean=mean(height), stError=sd(height)/sqrt(length(height)))
meanHeights <- transform(meanHeights, lower=mean-stError, upper=mean+stError)

print(
  ggplot(data=meanHeights) +
    geom_bar(mapping=aes(x=genotype, y=mean), fill="coral", stat="identity") +
    geom_errorbar(mapping=aes(x=genotype, ymax=upper, ymin=lower), width=0.5)
)

# To separate height data by genotype and group ---------------------------------------------------
meanHeightsGrp <- ddply(heightData, c("genotype", "group"), summarise, mean=mean(height), stError=sd(height)/sqrt(length(height)))
meanHeightsGrp <- transform(meanHeightsGrp, lower=mean-stError, upper=mean+stError)

print(
  ggplot(data=meanHeightsGrp) +
    geom_bar(mapping=aes(x=genotype, y=mean, fill=group), position="dodge", stat="identity") +
    geom_errorbar(mapping=aes(x=genotype, ymax=upper, ymin=lower, group=group), position=position_dodge(0.9), width=0.5)
)