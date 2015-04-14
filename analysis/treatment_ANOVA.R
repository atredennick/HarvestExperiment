##  Script to run ANOVAs for treatment effects of fire
##    and herbivory on tree regrowth after harvest.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.6.2015

##  This code is provided 'as is.' Any attempt to use this
##    code for different data is done so at the users own risk.

######
###### TODO: Check that mortality carries through based on NA regrowth for
######        non-NA initial biomass.
######

# Clear everything
rm(list=ls()) 

####
####  Load libraries --------------------------------------------------------
####
#install.packages(c('ggplot2','plyr','ggthemes','lme4'))
library('ggplot2')
library('plyr')
library('ggthemes')
library('lme4')


####
#### Read in data files ------------------------------------------------------
####
# Get treatment data and calculate mortality proportions
treatment_data <- read.csv("../data/tree_regrowth.csv")
#Loop over site and years
get_mort <- function(dataframe){
  num_dead <- length(which(is.na(dataframe$wood_wet_weight_kg)==TRUE))
  num_total <- nrow(dataframe)
  prop_mort <- num_dead / num_total
  return(c(prop_mort, num_dead, num_total))
}
TienAll <- get_mort(subset(treatment_data, site=="Tiendega"))
Tien2011 <- get_mort(subset(treatment_data, site=="Tiendega" & harvest_year==2011))
Tien2013 <- get_mort(subset(treatment_data, site=="Tiendega" & harvest_year==2013))
LakAll <- get_mort(subset(treatment_data, site=="Lakamane"))
Lak2011 <- get_mort(subset(treatment_data, site=="Lakamane" & harvest_year==2011))
Lak2013 <- get_mort(subset(treatment_data, site=="Lakamane" & harvest_year==2013))
mort_data <- data.frame(site = rep(c("Tiendega", "Lakamane"), each=3),
                        year = rep(c("2011", "2013", "average"), times=2),
                        prop_mort = c(Tien2011[1], Tien2013[1], TienAll[1],
                                      Lak2011[1], Lak2013[1], LakAll[1]),
                        num_dead = c(Tien2011[2], Tien2013[2], TienAll[2],
                                     Lak2011[2], Lak2013[2], LakAll[2]),
                        num_total = c(Tien2011[3], Tien2013[3], TienAll[3],
                                      Lak2011[3], Lak2013[3], LakAll[3]))
write.csv(mort_data, "../results/proportion_mortality.csv")

# Convert treatment data to dry weights
wet_dry_data <- read.csv("../data/wet_dry_subsamples.csv")
wet_dry_data$ratio <- with(wet_dry_data, dry_wood_weight_g / wet_wood_weight_g)
mean_wet_dry_species <- ddply(wet_dry_data, .(species_name, year), summarise,
                              mean_wet_dry = round(mean(ratio),2))

#set Lakamane 2013 wet:dry to 0.5 (see supplementary text)
mean_wet_dry_species <- rbind(mean_wet_dry_species,
                              c("Combretum glutinosum", 2013, 0.5))

treatment_data <- merge(treatment_data, mean_wet_dry_species, 
                        by.x=c("species_name", "harvest_year"),
                        by.y=c("species_name", "year"))
treatment_data$wood_dry_weight_g <- with(treatment_data,
                                         (wood_wet_weight_kg*as.numeric(mean_wet_dry))*1000)
allD <- treatment_data
# allD$harvest_year <- as.character(allD$harvest_year)
plot_info_data <- read.csv("../data/ssde_harvest_trees.csv")
plot_info_data$plot_id <- substring(plot_info_data$Plot, 4, 4)


####
#### Set up data for analysis -------------------------------------
####
initial_biomass <- read.csv("../data/tree_initial_biomass.csv")
initial_biomass <- initial_biomass[,c("site", "treatment_code", "tree_id", "initial_biomass_g")]
allD <- merge(allD, initial_biomass, by=c("site", "treatment_code", "tree_id"))
allD$harvest_year <- as.numeric(allD$harvest_year)
allD$rrg <- with(allD, 
                 (log(wood_dry_weight_g + initial_biomass_g) - log(initial_biomass_g)) / (harvest_year-2010))
#add in f/h factors
fire <- c("y", "y", "n", "n")
herb <- c("y", "n", "y", "n")
allD$fire <- NA
allD$herb <- NA
for(i in 1:length(fire)){
  tmp <- which(allD$treatment_code==i)
  allD[tmp,"fire"] <- fire[i]
  allD[tmp,"herb"] <- herb[i]
}

#Remove NAs
allDnoNA <- subset(allD, is.na(rrg)==FALSE)

# add in plot info
# get rid of extra treatment columns
drop <- "Treatment"
plot_info_data <- plot_info_data[,!(names(plot_info_data) %in% drop)]
test <- merge(allDnoNA,plot_info_data, by.x=c("site","treatment_code","tree_id"), by.y=c("Site","Treatment.Code","TreeID"))
# make sure the merge worked correctly
test$Year.Regrowth.Cut <- test$Year.Regrowth.Cut+2010
test$Year.Regrowth.Cut[which(test$Year.Regrowth.Cut==2012)] <- 2013
sum(with(test, Year.Regrowth.Cut-harvest_year)) #should be 0
# remove extra columns
drops <- c("Site.Code","Plot","Year.Regrowth.Cut","Species")
test <- test[,!(names(test) %in% drops)]
names(test)[which(names(test)=="species_name")] <- "Species"
allDnoNA <- test
# log transform to meet normality and homoscedasticity assumptions
allDnoNA$rrg <- log(allDnoNA$rrg)


####
#### Treatment ANOVAs ---------------------------------------------
####
# set up storage matrix for results
facANOVA <- matrix(nrow=4, ncol=6) 

#Tiendega ANOVAs
facANOVA[1, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Tiendega" & harvest_year==2011)))[1:3,4], 3)
facANOVA[1, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Tiendega" & harvest_year==2011)))[1:3,5], 3)
facANOVA[2, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Tiendega" & harvest_year==2013)))[1:3,4], 3)
facANOVA[2, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Tiendega" & harvest_year==2013)))[1:3,5], 3)
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & harvest_year==2011)))
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Tiendega" & Year.Harvested==2013)))

#Lakamane ANOVAs

facANOVA[3, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Lakamane" & harvest_year==2011)))[1:3,4], 3)
facANOVA[3, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Lakamane" & harvest_year==2011)))[1:3,5], 3)
facANOVA[4, c(1,3,5)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Lakamane" & harvest_year==2013)))[1:3,4], 3)
facANOVA[4, c(2,4,6)] <- round(anova(lm(rrg~fire*herb, data=subset(allDnoNA, site=="Lakamane" & harvest_year==2013)))[1:3,5], 3)
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2011)))
# plot(lm(rrg~fire*herb, data=subset(allDnoNA, Site=="Lakamane" & Year.Harvested==2013)))

#print results to check
facANOVA

#format results for output
facANOVAD <- as.data.frame(facANOVA)
colnames(facANOVAD) <- c("F_fire", "P_fire", "F_herb", "P_herb", "F_interact", "P_interact")
rownames(facANOVAD) <- c("Tiendega_2011", "Tiendega_2013", 
                        "Lakamane_2011", "Lakamane_2013")

#write results file (Table 1)
write.csv(facANOVAD, "../results/factorial_ANOVA_results.csv")


####
#### Some extra tests -----------------------------------------
####
#Test for year x treament interactions
anova(lm(rrg~fire*herb*as.factor(harvest_year), data=subset(allDnoNA, site=="Tiendega")))
anova(lm(rrg~fire*herb*as.factor(harvest_year), data=subset(allDnoNA, site=="Lakamane")))

#Test for plot effects
#Tiendega
mix1 <- lmer(rrg ~ fire*herb + (1|plot_id), data=subset(allDnoNA, site=="Tiendega" & harvest_year==2011))
nomix1 <- lm(rrg~fire*herb, data=subset(allDnoNA, site=="Tiendega" & harvest_year==2011))
AIC(mix1, nomix1)
mix1 <- lmer(rrg ~ fire*herb + (1|plot_id), data=subset(allDnoNA, site=="Tiendega" & harvest_year==2013))
nomix1 <- lm(rrg~fire*herb, data=subset(allDnoNA, site=="Tiendega" & harvest_year==2013))
AIC(mix1, nomix1)
#Lakamane
mix1 <- lmer(rrg ~ fire*herb + (1|plot_id), data=subset(allDnoNA, site=="Lakamane" & harvest_year==2011))
nomix1 <- lm(rrg~fire*herb, data=subset(allDnoNA, site=="Lakamane" & harvest_year==2011))
AIC(mix1, nomix1)
mix1 <- lmer(rrg ~ fire*herb + (1|plot_id), data=subset(allDnoNA, site=="Lakamane" & harvest_year==2013))
nomix1 <- lm(rrg~fire*herb, data=subset(allDnoNA, site=="Lakamane" & harvest_year==2013))
AIC(mix1, nomix1)


####
#### Pooled treatments ANOVAs -----------------------------------
####
#set up storage matrix for results
pooledANOVA <- matrix(nrow=2, ncol=4)
#by year
pooledANOVA[,1] <- round(as.numeric(anova(aov(rrg~site, data=subset(allDnoNA, harvest_year==2011)))[1,4:5]),3)
TukeyHSD(aov(rrg~site, data=subset(allDnoNA, harvest_year==2011)))
pooledANOVA[,2] <- round(as.numeric(anova(aov(rrg~site, data=subset(allDnoNA, harvest_year==2013)))[1,4:5]),3)
TukeyHSD(aov(rrg~site, data=subset(allDnoNA, harvest_year==2013)))

#by site
pooledANOVA[,3] <- round(as.numeric(anova(aov(rrg~as.factor(harvest_year), data=subset(allDnoNA, site=="Tiendega")))[1,4:5]),3)
TukeyHSD(aov(rrg~as.factor(harvest_year), data=subset(allDnoNA, site=="Tiendega")))
pooledANOVA[,4] <- round(as.numeric(anova(aov(rrg~as.factor(harvest_year), data=subset(allDnoNA, site=="Lakamane")))[1,4:5]),3)
TukeyHSD(aov(rrg~as.factor(harvest_year), data=subset(allDnoNA, site=="Lakamane")))

#format results for output
poolD <- as.data.frame(pooledANOVA)
colnames(poolD) <- c("SiteDiff2011", "SiteDiff2013", "YearDiffTiendega", "YearDiffLakamane")
rownames(poolD) <- c("F", "P")

#write results file
write.csv(poolD, "../results/pooled_treatment_results.csv")



####
#### Make figure 3 --------------------------------------------------
####
#change rrg back to non-log version
allDnoNA <- subset(allD, is.na(rrg)==FALSE)

#Figure 3: Factorial ANOVA
pdf("../results/Figure3.pdf", height=5, width=7)
print(
ggplot(data=allDnoNA, aes(x=fire, y=rrg, fill=herb))+
  geom_boxplot(outlier.shape = 1)+
  facet_grid(harvest_year~site)+
  ylab(expression(paste("Relative Regrowth Rate (g ", g^-1, " ", y^-1,")")))+
  xlab("Fire Treatment")+
  scale_x_discrete(labels=c("No Fire", "Fire"))+
  scale_fill_manual(labels=c("No Large Hebivores", "Large Herbivores"), values=c("grey40", "grey80"),
                    name="Herbivore Treatment")+
  theme_bw()
)
dev.off()

# plot by site only for each year
# ggplot(data=allDnoNA, aes(x=Site, y=rrg))+
#   geom_boxplot(outlier.shape = 1, fill="grey", width=0.5)+
#   ylab(expression(paste("Relative Regrowth Rate (g ", g^-1, " ", y^-1,")")))+
#   facet_grid(Year.Harvested~.)+
#   theme_bw()

