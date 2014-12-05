rm(list=ls(all=TRUE)) #clear everything, just to be safe

##TO DO: Estimate 2013 biomass change for fh adult trees; compare to 2013 regrowth. Then DONE!

#load some libraries
#if you need to install these packages, remove the '#' from the next line and run them
#install.packages(c("plyr", "ggplot2"))
library(plyr)
library(ggplot2)
library(ggthemes)
library(lme4)



wdVec <- seq(0.1, 1, 0.01)
p1 <- numeric(length(wdVec))
p2 <- numeric(length(wdVec))

for(j in 1:length(wdVec)){
  y3D <- read.csv("Regrowth_2013.csv")
  y3D <- subset(y3D, Site!="Tiorola") #remove Tiorola
  #wet:dry ratios for 2013 (Lakamane not estimated); Tiendega from 2013 harvest subsamples
  wdRatios <- data.frame(Site = c("Tiendega", "Lakamane"))
  wdRatios$Ratio <- c(wdVec[j], 0.44)
  tmpD <- merge(y3D[,1:7], wdRatios, all.x=TRUE)
  regrow2013_kg <- tmpD[,6]*tmpD[,8]
#   rm(y3D, wdRatios, tmpD)
  
  ##Bring in all data
  allD <- read.csv("Regrowth_2011-2013.csv")
  allD <- subset(allD, Site!="Tiorola")
  tmp <- which(allD$Year.Harvested==2013)
  allD[tmp, 6] <- regrow2013_kg
  rm(tmp, regrow2013_kg)
  
  #Rescale things
  allD$Initial.Biomass..kg. <- allD$Initial.Biomass..g./1000
  allD$rrg <- (log(allD$TOTAL.BIOMASS..kg.+allD$Initial.Biomass..kg.)-log(allD$Initial.Biomass..kg.))/(allD$Year.Harvested-2010)
  allD$logrrg <- log(allD$rrg)
  
  #add in f/h factors
  fire <- c("y", "y", "n", "n")
  herb <- c("y", "n", "y", "n")
  allD$fire <- NA
  allD$herb <- NA
  for(i in 1:length(fire)){
    tmp <- which(allD$Treatment==i)
    allD[tmp,"fire"] <- fire[i]
    allD[tmp,"herb"] <- herb[i]
  }
  
  #Remove NAs
  allDnoNA <- subset(allD, is.na(rrg)==FALSE)
  
  # #Pooled treatment, site effects by year
  mod1 <- aov(logrrg~as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Lakamane"))
  quickSub <- subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2013)
  print(mean(quickSub$rrg))
  anova(mod1)
  TukeyHSD(mod1)
  
  mod2 <- aov(logrrg~Site, data=subset(allDnoNA, Year.Harvested==2013))
  anova(mod2)
  TukeyHSD(mod2)
  
  p1[j] <- summary(mod1)[[1]][["Pr(>F)"]][[1]]
  p2[j] <- summary(mod2)[[1]][["Pr(>F)"]][[1]]
}

# par(family="serif")
pdf(file = "../Pvalue_TestWetDry.pdf", width = 7, height = 6)
plot(wdVec, p1, xlab="dry:wet ratio", ylab="p-value", pch=21, col="black", bg="grey", xlim=c(0,1))
points(wdVec, p2, pch=21, col="black", bg="white")
abline(h=0.1, lwd=4)
points(0.8, 0.11, pch=22, col="white", bg="white", cex=10)
text(0.8, 0.1, "p = 0.1")
abline(v=0.5, lty=2, lwd=4)
# points(0.5, 0.8, pch=22, col="white", bg="white", cex=7)
# text(0.4, 0.9, "dry:wet = 0.5")
legend(0, 1, legend = c("Year diff. (in Lakamane)", "Site diff. (in 2013)"), pt.bg = c("grey", "white"), pch=21, bty="n")
dev.off()



