# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Script to reproduce statistical results and figures for Tredennick et al. #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Contact: Andrew Tredennick 
# Email: atredenn@gmail.com

rm(list=ls(all=TRUE)) #clear everything, just to be safe



####
#### Load libraries
####
#if you need to install these packages, remove the '#' from the next line and run them (need internet connection)
#install.packages(c("plyr", "ggplot2", "ggthemes", "lme4"))
library(plyr)
library(ggplot2)
library(ggthemes)
library(lme4)



####
#### Read in data file
####
allD <- read.csv("AllData_ForAnalysis_2011-2013.csv", na.strings = "NA")



####
#### Get percent mortality
####
Tien1 <- length(which(is.na(subset(allD, Site=="Tiendega" & Year.Harvested=="2011"))[,7]==TRUE))/nrow(subset(allD, Site=="Tiendega" & Year.Harvested=="2011"))
Tien2 <- length(which(is.na(subset(allD, Site=="Tiendega" & Year.Harvested=="2013"))[,7]==TRUE))/nrow(subset(allD, Site=="Tiendega" & Year.Harvested=="2013"))
Lak1 <- length(which(is.na(subset(allD, Site=="Lakamane" & Year.Harvested=="2011"))[,7]==TRUE))/nrow(subset(allD, Site=="Lakamane" & Year.Harvested=="2011"))
Lak2 <- length(which(is.na(subset(allD, Site=="Lakamane" & Year.Harvested=="2013"))[,7]==TRUE))/nrow(subset(allD, Site=="Lakamane" & Year.Harvested=="2013"))
TienTot <- length(which(is.na(subset(allD, Site=="Tiendega"))[,7]==TRUE))/nrow(subset(allD, Site=="Tiendega"))
LakTot <- length(which(is.na(subset(allD, Site=="Lakamane"))[,7]==TRUE))/nrow(subset(allD, Site=="Lakamane"))

#print results as percentages
Tien1*100
Tien2*100
Lak1*100
Lak2*100
TienTot*100
LakTot*100



####
#### Set up data for analysis
####
#Remove NAs
allDnoNA <- subset(allD, is.na(rrg)==FALSE)
#log transform to meet normality and homoscedasticity assumptions
allDnoNA$rrg <- log(allDnoNA$rrg)



####
#### Treatment ANOVAs
####
#set up storage matrix for results
facANOVA <- matrix(nrow=6, ncol=6) 

#Tiendega ANOVAs
facANOVA[1, c(1,3,5)] <- round(anova(lmer(rrg ~ fire*herb + (1|Site), data=subset(allDnoNA, Year.Harvested==2011)))[1:3,4],3)
facANOVA[2, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2011)))[1:3,4], 3)
facANOVA[2, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2011)))[1:3,5], 3)
facANOVA[3, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2013)))[1:3,4], 3)
facANOVA[3, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2013)))[1:3,5], 3)
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2011)))
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2013)))

#Lakamane ANOVAs
facANOVA[4, c(1,3,5)] <- round(anova(lmer(rrg ~ fire*herb + (1|Site), data=subset(allDnoNA, Year.Harvested==2011)))[1:3,4],3)
facANOVA[5, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2011)))[1:3,4], 3)
facANOVA[5, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2011)))[1:3,5], 3)
facANOVA[6, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2013)))[1:3,4], 3)
facANOVA[6, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2013)))[1:3,5], 3)
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2011)))
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2013)))

#print results to check
facANOVA

#format results for output
facANOVAD <- as.data.frame(facANOVA)
colnames(facANOVAD) <- c("F_fire", "P_fire", "F_herb", "P_herb", "F_interact", "P_interact")
rownames(facANOVAD) <- c("Tiendega_Mixed", "Tiendega_2011", "Tiendega_2013", 
                        "Lakamane_Mixed", "Lakamane_2011", "Lakamane_2013")

#write results file
write.csv(facANOVAD, "factorialANOVA_results.csv")

#Test for year x treament interactions
anova(lm(rrg~fire*herb*as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Tiendega")))
anova(lm(rrg~fire*herb*as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Lakamane")))



####
#### Pooled treatments ANOVAs
####
#set up storage matrix for results
pooledANOVA <- matrix(nrow=2, ncol=4)
#by year
pooledANOVA[,1] <- round(as.numeric(anova(aov(rrg~Site, data=subset(allDnoNA, Year.Harvested==2011)))[1,4:5]),3)
TukeyHSD(aov(rrg~Site, data=subset(allDnoNA, Year.Harvested==2011)))
pooledANOVA[,2] <- round(as.numeric(anova(aov(rrg~Site, data=subset(allDnoNA, Year.Harvested==2013)))[1,4:5]),3)
TukeyHSD(aov(rrg~Site, data=subset(allDnoNA, Year.Harvested==2013)))

#by site
pooledANOVA[,3] <- round(as.numeric(anova(aov(rrg~as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Tiendega")))[1,4:5]),3)
TukeyHSD(aov(rrg~as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Tiendega")))
pooledANOVA[,4] <- round(as.numeric(anova(aov(rrg~as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Lakamane")))[1,4:5]),3)
TukeyHSD(aov(rrg~as.factor(Year.Harvested), data=subset(allDnoNA, Site=="Lakamane")))

#format results for output
poolD <- as.data.frame(pooledANOVA)
colnames(poolD) <- c("SiteDiff2011", "SiteDiff2013", "YearDiffTiendega", "YearDiffLakamane")
rownames(poolD) <- c("F", "P")

#write results file
write.csv(poolD, "pooledTrtResults.csv")



####
#### Make figures 3 and 4
####
#change rrg back to non-log version
allDnoNA <- subset(allD, is.na(rrg)==FALSE)

#Figure 3: Factorial ANOVA
ggplot(data=allDnoNA, aes(x=fire, y=rrg, fill=herb))+
  geom_boxplot(outlier.shape = 1)+
  facet_grid(Year.Harvested~Site)+
  ylab(expression(paste("Relative Regrowth Rate (g ", g^-1, " ", y^-1,")")))+
  xlab("Fire Treatment")+
  scale_x_discrete(labels=c("No Fire", "Fire"))+
  scale_fill_manual(labels=c("No Large Hebivores", "Large Herbivores"), values=c("grey40", "grey80"),
                    name="Herbivore Treatment")+
  theme_bw()

# Figure 4: plot by site only for each year
ggplot(data=allDnoNA, aes(x=Site, y=rrg))+
  geom_boxplot(outlier.shape = 1, fill="grey", width=0.5)+
  ylab(expression(paste("Relative Regrowth Rate (g ", g^-1, " ", y^-1,")")))+
  facet_grid(Year.Harvested~.)+
  theme_bw()

