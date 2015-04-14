##  Script to analyze height response after harvest. We use simple
##    ANOVAs to test for treatment effects and then pool the data
##    to look at the probability of reaching escape heights based on
##    initial tree biomass.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.7.2015


#clear everything
rm(list=ls()) 


####
####  Load libraries ---------------------------------------------
####
library('ggplot2')
library('nlme')
library('plyr')
library('reshape2')
library('ggthemes')
library('gridExtra')


####
####  Read in 2013 data and plot by treatment and site -----------
####
allD <- read.csv("../data/heights2013.csv")

# Nice boxplots of heights by treatment and site
g1 <- ggplot(allD, aes(x=as.factor(Treatment), y=Height_m))+
  geom_hline(aes(yintercept=1), linetype=3)+
  geom_hline(aes(yintercept=2), linetype=2)+
  geom_boxplot(fill="grey80", outlier.shape=19, alpha=0.9, na.rm=TRUE)+
  facet_grid(Site~.)+
  ylab("Height after three years (m)")+
  xlab("Treatment")+
  scale_x_discrete(labels=c("FH", "Fh", "fH", "fh"))+
  ggtitle("A")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))

# Quick ANOVA to test for treatment effects
#Tiendega
anova(lm(log(Height_m)~as.factor(Treatment), data=subset(allD, Site=="Tiendega")))
mod <- aov(log(Height_m)~as.factor(Treatment), data=subset(allD, Site=="Tiendega"))
anova(mod)
TukeyHSD(mod)

#Lakamane
anova(lm(log(Height_m)~as.factor(Treatment), data=subset(allD, Site=="Lakamane")))
mod <- aov(log(Height_m)~as.factor(Treatment), data=subset(allD, Site=="Lakamane"))
anova(mod)
TukeyHSD(mod)

#Test for site differences -- only N of 2, so ignore this...
# site_diffs <- anova(lm(log(Height_m)~Site, data=allD))
# site_height_diffs <- data.frame(F_value=site_diffs[[4]][1],
#                                 P_value=site_diffs[[5]][1])

# Save results
site_height_stats <- ddply(allD, .(Site), summarise,
                           avg_height = mean(Height_m, na.rm = TRUE),
                           sd_height = sd(Height_m, na.rm = TRUE))
write.csv(site_height_stats, "../results/site_height_stats.csv")


####
#### Logistic models for height ---------------------------------------------------
####

###
### TIENDEGA
###

##  2 meter height threshold
hts <- data.frame(Height = allD[!is.na(allD$Height_m),"Height_m"],
                  Site = allD[!is.na(allD$Height_m),"Site"],
                  Treatment = allD[!is.na(allD$Height_m),"Treatment"],
                  Initial = allD[!is.na(allD$Height_m),"InitBiomass_g"])
htThresh <- 2
hts$Above <- NA
hts[hts$Height >= htThresh,ncol(hts)] <- "Above"
hts[hts$Height < htThresh,ncol(hts)] <- "Below"
hts$Rec <- 0
hts[which(hts$Above == "Above"), ncol(hts)] <- 1
hts <- hts[hts$Initial < 250000,]

rec.out2<-glm(Rec~Initial,data=subset(hts, Site=="Tiendega"),family="binomial")
tiendega_2meter_glm <- summary(rec.out2)
newdata2 <- with(hts, data.frame(Initial = seq(0,250000, by=1000)))
predS2 <- as.data.frame(predict(rec.out2, newdata=newdata2, type="link", se=TRUE))
predS2 <- within(predS2, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
predS2$Initial <- seq(0,250000, by=1000)


##  3 meter height threshold
hts <- data.frame(Height = allD[!is.na(allD$Height_m),"Height_m"],
                  Site = allD[!is.na(allD$Height_m),"Site"],
                  Treatment = allD[!is.na(allD$Height_m),"Treatment"],
                  Initial = allD[!is.na(allD$Height_m),"InitBiomass_g"])
htThresh <- 3
hts$Above <- NA
hts[hts$Height >= htThresh,ncol(hts)] <- "Above"
hts[hts$Height < htThresh,ncol(hts)] <- "Below"
hts$Rec <- 0
hts[which(hts$Above == "Above"), ncol(hts)] <- 1
hts <- hts[hts$Initial < 250000,]

rec.out2<-glm(Rec~Initial,data=subset(hts, Site=="Tiendega"),family="binomial")
tiendega_3meter_glm <- summary(rec.out2)
newdata2 <- with(hts, data.frame(Initial = seq(0,250000, by=1000)))
predS3 <- as.data.frame(predict(rec.out2, newdata=newdata2, type="link", se=TRUE))
predS3 <- within(predS3, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
predS3$Initial <- seq(0,250000, by=1000)

g2 <- ggplot()+
#   geom_ribbon(data=predS1, aes(ymin=0, ymax=PredictedProb, x=Initial/100), fill="grey90")+
#   geom_line(data=predS1, aes(x = Initial/100, y = PredictedProb), color="black", size=1, linetype=1)+
  geom_ribbon(data=predS2, aes(ymin=0, ymax=PredictedProb, x=Initial/100), fill="grey70", color="white")+
  geom_line(data=predS2, aes(x = Initial/100, y = PredictedProb), color="black", size=1, linetype=2)+
  geom_ribbon(data=predS3, aes(ymin=0, ymax=PredictedProb, x=Initial/100), fill="white", color="white")+
  geom_line(data=predS3, aes(x = Initial/100, y = PredictedProb), color="black", size=1, linetype=3)+
  xlab("Initial Biomass (kg)") + ylab("Probability of Achieving Escape Height")+
  geom_text(aes(x=2000, y=0.1), label = "Tiendega", size=6)+
  ggtitle("B")+
  theme_few()+
  theme(plot.title = element_text(hjust = 1))


###
### Lakamane
###
##  2 meter threshold
hts <- data.frame(Height = allD[!is.na(allD$Height_m),"Height_m"],
                  Site = allD[!is.na(allD$Height_m),"Site"],
                  Treatment = allD[!is.na(allD$Height_m),"Treatment"],
                  Initial = allD[!is.na(allD$Height_m),"InitBiomass_g"])
htThresh <- 2
hts$Above <- NA
hts[hts$Height >= htThresh,ncol(hts)] <- "Above"
hts[hts$Height < htThresh,ncol(hts)] <- "Below"
hts$Rec <- 0
hts[which(hts$Above == "Above"), ncol(hts)] <- 1
hts <- hts[hts$Initial < 250000,]

rec.out2<-glm(Rec~Initial,data=subset(hts, Site=="Lakamane"),family="binomial")
lakamane_2meter_glm <- summary(rec.out2)
newdata2 <- with(hts, data.frame(Initial = seq(0,250000, by=1000)))
predS2 <- as.data.frame(predict(rec.out2, newdata=newdata2, type="link", se=TRUE))
predS2 <- within(predS2, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
predS2$Initial <- seq(0,250000, by=1000)


##  3 meter threshold
hts <- data.frame(Height = allD[!is.na(allD$Height_m),"Height_m"],
                  Site = allD[!is.na(allD$Height_m),"Site"],
                  Treatment = allD[!is.na(allD$Height_m),"Treatment"],
                  Initial = allD[!is.na(allD$Height_m),"InitBiomass_g"])
htThresh <- 3
hts$Above <- NA
hts[hts$Height >= htThresh,ncol(hts)] <- "Above"
hts[hts$Height < htThresh,ncol(hts)] <- "Below"
hts$Rec <- 0
hts[which(hts$Above == "Above"), ncol(hts)] <- 1
hts <- hts[hts$Initial < 250000,]

rec.out2<-glm(Rec~Initial,data=subset(hts, Site=="Lakamane"),family="binomial")
lakamane_3meter_glm <- summary(rec.out2)
newdata2 <- with(hts, data.frame(Initial = seq(0,250000, by=1000)))
predS3 <- as.data.frame(predict(rec.out2, newdata=newdata2, type="link", se=TRUE))
predS3 <- within(predS3, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
predS3$Initial <- seq(0,250000, by=1000)

g3 <- ggplot()+
#   geom_ribbon(data=predS1, aes(ymin=0, ymax=PredictedProb, x=Initial/100), fill="grey90")+
#   geom_line(data=predS1, aes(x = Initial/100, y = PredictedProb), color="black", size=1, linetype=1)+
  geom_ribbon(data=predS2, aes(ymin=0, ymax=PredictedProb, x=Initial/100), fill="grey70", color="white")+
  geom_line(data=predS2, aes(x = Initial/100, y = PredictedProb), color="black", size=1, linetype=2)+
  geom_ribbon(data=predS3, aes(ymin=0, ymax=PredictedProb, x=Initial/100), fill="white", color="white")+
  geom_line(data=predS3, aes(x = Initial/100, y = PredictedProb), color="black", size=1, linetype=3)+
  xlab("Initial Biomass (kg)") + ylab("Probability of Achieving Escape Height")+
  geom_text(aes(x=2000, y=0.1), label = "Lakamane", size=6)+
#   geom_text(aes(x=150, y=0.8), label = "1 m", size=4)+
  geom_text(aes(x=500, y=0.5), label = "2 m", size=4)+
  geom_text(aes(x=1000, y=0.2), label = "3 m", size=4)+
  ggtitle("C")+
  theme_few()+
  theme(plot.title = element_text(hjust = 1))


####
####  Make Figure 4 -----------------------------------------------------------------------
####
pdf("../results/Figure4.pdf", width = 4, height = 10)
grid.arrange(g1, g2, g3, nrow=3, ncol=1)
dev.off()


####
####  Save GLM output ----------------------------------------------------------------------
####
row.names(tiendega_2meter_glm$coefficients) <- NULL
row.names(tiendega_3meter_glm$coefficients) <- NULL
row.names(lakamane_2meter_glm$coefficients) <- NULL
row.names(lakamane_3meter_glm$coefficients) <- NULL
out_glms <- rbind(tiendega_2meter_glm$coefficients,
                  tiendega_3meter_glm$coefficients,
                  lakamane_2meter_glm$coefficients,
                  lakamane_3meter_glm$coefficients)
out_glms <- as.data.frame(out_glms)
out_glms$site <- rep(c("Tiendega","Lakamane"), each=4)
out_glms$height_threshold <- rep(c("2","2","3","3"), times=2)
write.csv(out_glms, "../results/glm_height_threshold_results.csv")

