##  This script runs statistical analyses that involve comparing
##    2013 Lakamane data to non-2013 Lakamane data since the wet:dry
##    ratio data for Lakamane in 2013 was lost in a storm. The idea
##    is to use a range of wet:dry ratios, apply those to the wet biomass
##    data, run the statistical analyses, and then see how sensitive the
##    P-values are to the wet:dry ratio. If the tests are insensitive to
##    to the ratio, we can be sure of our conclusions. If not, our conclusions
##    are suspect, or must be taken with a grain of salt. Importantly,
##    the main results in the paper on treatment effects by year and site
##    are impervious to the wet:dry ratio since it is averaged by site and
##    year.

##  There are two tests that could be affected by biased/wrong wet:dry ratio in
##    in 2013:
##          (1) Test between control tree growth and 2013 protected regrowth
##          (2) Test between 2011 regrowth and 2013 regrowth pooled over treatments

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.7.2015


rm(list=ls(all=TRUE)) #clear everything, just to be safe
library('plyr')
library('ggplot2')

####
####  Set up wet:dry vector and bring in data -----------------------------------
####
wdVec <- seq(0.1, 1, 0.01)



#Loop over wet:dry ratios and perform statistical tests
p_harvcontrol <- numeric(length(wdVec))
p_lakyear <- numeric(length(wdVec))
for(wd_now in 1:length(wdVec)){
  # Bring in treatment regrowth data and wet/dry weights
  treatment_data <- read.csv("../data/tree_regrowth.csv")
  wet_dry_data <- read.csv("../data/wet_dry_subsamples.csv")
  wet_dry_data$ratio <- with(wet_dry_data, dry_wood_weight_g / wet_wood_weight_g)
  mean_wet_dry_species_main <- ddply(wet_dry_data, .(species_name, year), summarise,
                                     mean_wet_dry = round(mean(ratio),2))
  
  # Bring in initial biomass data
  initial_biomass <- read.csv("../data/tree_initial_biomass.csv")
  initial_biomass <- initial_biomass[,c("site", "treatment_code", 
                                        "tree_id", "initial_biomass_g")]
  
  # Bring in control tree growth data (results from allometry predictions)
  control_data <- read.csv("../results/control_predicted_biomass.csv")
  control_data <- control_data[,c("Site", "Species", "Year", "rgr")]
  control_data$Treatment <- "C"
  set_to_na <- which(control_data$rgr==0)
  control_data$rgr[set_to_na] <- NA
  #set Lakamane 2013 wet:dry to desired value from vector
  mean_wet_dry_species <- rbind(mean_wet_dry_species_main,
                                c("Combretum glutinosum", 2013, wdVec[wd_now]))
  treatment_data <- merge(treatment_data, mean_wet_dry_species, 
                          by.x=c("species_name", "harvest_year"),
                          by.y=c("species_name", "year"))
  treatment_data$wood_dry_weight_g <- with(treatment_data,
                                           (wood_wet_weight_kg*as.numeric(mean_wet_dry))*1000)
  allD <- treatment_data
  allD <- merge(allD, initial_biomass, by=c("site", "treatment_code", "tree_id"))
  allD$harvest_year <- as.numeric(allD$harvest_year)
  allD$rrg <- with(allD, 
                   (log(wood_dry_weight_g + initial_biomass_g) - log(initial_biomass_g)) / (harvest_year-2010))
  
  # Subset fh data
  fh_data <- subset(allD, treatment_code==4)
  fh_data$treatment <- "H"
  fh_data <- fh_data[,c("site", "species_name", "harvest_year", "rrg", "treatment")]
  
  # Combine the two datasets
  colnames(fh_data) <- colnames(control_data)
  allD <- rbind(control_data, fh_data)
  colnames(allD)[which(colnames(allD)=="rgr")] <- "rrg"
  
  # Transform the regrowth
  allD$rrgNorm <- allD$rrg
  allD$rrg <- log(allD$rrg) # log transform to meet normality
  # testing of normality assumptions not shown
  allDnoNA <- subset(allD, is.na(rrg)==FALSE)
  mod1 <- aov(rrg~as.factor(Treatment)*as.factor(Year), data=subset(allDnoNA, Site=="Lakamane"))
  p_harvcontrol[wd_now] <- summary(mod1)[[1]][["Pr(>F)"]][[1]]
  
  
  # Now do year effect for Lakamane with pooled treatment data
  allD <- treatment_data
  allD <- merge(allD, initial_biomass, by=c("site", "treatment_code", "tree_id"))
  allD$harvest_year <- as.numeric(allD$harvest_year)
  allD$rrg <- with(allD, 
                   (log(wood_dry_weight_g + initial_biomass_g) - log(initial_biomass_g)) / (harvest_year-2010))
  allDnoNA <- subset(allD, is.na(rrg)==FALSE)
  mod2 <- aov(log(rrg)~as.factor(harvest_year), data=subset(allDnoNA, site=="Lakamane"))
  p_lakyear[wd_now] <- summary(mod2)[[1]][["Pr(>F)"]][[1]]
}

####
####  Plot the results ------------------------------------------------
####
plot_data <- data.frame(ratio=rep(wdVec,2),
                        p_value=c(p_harvcontrol,p_lakyear),
                        test=rep(c("fh vs. Control", "2011 vs. 2013"), each=length(wdVec)))
pdf(file = "../results/FigureS2.pdf", width = 5, height = 3)
print(
  ggplot(plot_data, aes(x=ratio, y=p_value))+
  geom_point()+
  facet_wrap("test", scales="free")
)
dev.off()
