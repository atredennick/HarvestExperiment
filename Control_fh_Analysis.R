rm(list=ls()) 
library(ggplot2)
library(nlme)
library(plyr)
library(reshape2)

#Control vs. fh regrowth analysis
allD <- read.csv("controlGrowth.csv")
allD$rrgNorm <- allD$rrg
allD$rrg <- log(allD$rrg)

#Run some simple ANOVAs
mod1 <- aov(rrg~as.factor(Treatment)*Site, data=subset(allD, Year=="2011"))
anova(mod1)
TukeyHSD(mod1)

mod1 <- aov(rrg~as.factor(Treatment)*Site, data=subset(allD, Year=="2013"))
anova(mod1)
TukeyHSD(mod1)

mod1 <- aov(rrg~as.factor(Treatment)*as.factor(Year), data=subset(allD, Site=="Tiendega"))
anova(mod1)
T1 <- TukeyHSD(mod1)

mod2 <- aov(rrg~as.factor(Treatment)*as.factor(Year), data=subset(allD, Site=="Lakamane"))
anova(mod2)
T2 <- TukeyHSD(mod2)

#Effect size: average difference
avgC <- ddply(allD, .(Year, Site, Treatment), summarise,
                  avgRRGR = mean(rrgNorm, na.rm = TRUE))
lak1 <- subset(avgC, Year==2011 & Site=="Lakamane")[1,4]/subset(avgC, Year==2011 & Site=="Lakamane")[2,4]
lak2 <- subset(avgC, Year==2013 & Site=="Lakamane")[1,4]/subset(avgC, Year==2013 & Site=="Lakamane")[2,4]
tien1 <- subset(avgC, Year==2011 & Site=="Tiendega")[1,4]/subset(avgC, Year==2011 & Site=="Tiendega")[2,4]
tien2 <- subset(avgC, Year==2013 & Site=="Tiendega")[1,4]/subset(avgC, Year==2013 & Site=="Tiendega")[2,4]

mean(c(lak1, lak2))
mean(c(tien1, tien2))

#Make some plots
ggplot(allD, aes(x=as.factor(Year), y=rrgNorm, fill=Treatment))+
  geom_boxplot(outlier.shape = 1)+
  facet_grid(Site~.)+
  scale_fill_manual(labels=c("Control", "Harvested (fh)"), values=c("grey40", "grey80"),
                    name="")+
  ylab(expression(paste("Relative (Re)Growth Rate (g ", g^-1, " ", y^-1,")")))+
  xlab("Measurement Year")+
  theme_bw()
