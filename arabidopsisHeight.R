library(ggplot2)
library(plyr)

heightData <- read.csv(file="E:/Scutellaria/Col and A2 -1 combined.csv")
colnames(heightData)  <- c("genotype", "group", "height")
heightData[, 2] <- as.factor(heightData[, 2])

meanHeights <- ddply(heightData, c("genotype", "group"), summarise, mean=mean(height), stError=sd(height)/sqrt(length(height)))
meanHeights <- transform(meanHeights, lower=mean-stError, upper=mean+stError)

print(
ggplot(data=meanHeights) +
  geom_bar(mapping=aes(x=genotype, y=mean, fill=group), position="dodge", stat="identity") +
  geom_errorbar(mapping=aes(x=genotype, ymax=upper, ymin=lower, group=group), position=position_dodge(0.9), width=0.5)
)