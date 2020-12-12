library(ggplot2)
setwd("Documents/CIS4900/")

myvals <- read.csv("Select-real.txt", header = FALSE)

totalAvg <- mean(myvals[,2])
totalSe <- sd(myvals[,2])/sqrt(length(myvals[,2]))

data <- data.frame(
  name=myvals[,1],
  value=myvals[,2],
  se=myvals[,3],
  avg=totalAvg,
  avgSe=totalSe
)

ggplot(data) +
  geom_bar(aes(x=name, y=value), stat="identity", fill="grey", alpha=1) +
  geom_pointrange(aes(x=name, y=value, ymin=value-se, ymax=value+se), colour="black", alpha=0.9, size=0.4) +
  coord_cartesian(ylim=c(0.9,1)) +
  ggtitle("Average Proportion of Haplotypes Found for Values of N*\n in Real Data Example, Population Size 10,000") +
  xlab("Values of N*") +
  ylab("Proportion of Haplotypes Found") + 
  geom_abline(slope = 0, intercept = totalAvg, col="red")