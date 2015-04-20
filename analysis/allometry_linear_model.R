##  Script to test alternative allometric models,
##    choose best model based on AICc, and then predict
##    biomass given observed diameters. The model 
##    predicts log(biomass) in grams from diameter in 
##    centimeters.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.6.2015

##  This code is provided 'as is.' Any attempt to use this
##    code for different data is done so at the users own risk.

# Clear everything from the workspace
rm(list=ls())

####
####  Load libraries --------------------------------
####
#install.packages(c('ggplot2', 'nlme', 'plyr', 'reshape2')) # uncomment if you need these
library('ggplot2')
library('nlme')
library('plyr')
library('reshape2')
library('AICcmodavg')

####
####  Bring in data --------------------------------
####
allometry_data <- read.csv("../data/allometry_data_main_stems.csv")
names(allometry_data) <- tolower(names(allometry_data))

allometry_data <- allometry_data[allometry_data$aggregated_leaf_wt>-9000,]


Y <- allometry_data$aggregated_wood_wt
x <- allometry_data$diameter
species <- allometry_data$spp
tree <- allometry_data$tree


####
####  Compare models ---------------------------------
####
# Compare models with/without species random effect
#With species and tree random effect
mix <- lme(fixed=log(Y)~log(x), random = ~ 1 | species/tree, method="ML")
summary(mix)

# With species random effect
mix1 <- lme(fixed=log(Y)~log(x), random = ~ 1 | species, method="ML")
summary(mix1)

# Update to remove species effect
species_null <- (numeric(length(species)))
mix2 <- update(mix1, random = ~ 1 | species_null)
summary(mix2)

##Compare models
anova(mix, mix1)
anova(mix1, mix2)
AICc1 <- AICc(mix)
AICc2 <- AICc(mix1)
AICc3 <- AICc(mix2)
AICc1
AICc2
AICc3
#AICc scores are really close, use simpler model with no species random effect

# Save the AICc scores
aic_scores <- data.frame(model = c("spp-tree_mix", "spp_mix", "no_rand"),
                         aicc = c(AICc1, AICc2, AICc3))
write.csv(aic_scores, "../results/allometric_model_aics.csv")

####
####  Fit final model ------------------------------------
####
model <- lme(fixed=log(Y)~log(x), random = ~ 1 | species, method="ML")
summary(model)

newx <- rep(seq(min(x), max(x), by=0.1),times=2)
newspp <- rep(unique(allometry_data$spp), each=(length(newx)/2))
# Plot for supplementary material
pred_plot <- as.numeric(predict(model, newdata=data.frame(x=newx,species=newspp)))
pred_data <- data.frame(diameter=newx, biomass=pred_plot, species=newspp)
plot_data <- data.frame(diameter=x, biomass=Y,
                        species=allometry_data$spp)

pdf("../results/FigureS1.pdf", width = 6, height = 3.5)
print(
ggplot()+
  geom_point(data=plot_data, aes(log(diameter), log(biomass), color=species), size=3.5)+
  geom_point(data=plot_data, aes(log(diameter), log(biomass)),shape=1, size=3.5, color="black")+
  geom_line(data=pred_data, aes(log(diameter), biomass, color=species), size=1)+
  scale_color_manual(values=c("darkorange", "steelblue"))
)
dev.off()


####
####  Predict new values from fitted model ----------------
####
# Pull in data to predict initial tree biomass (from diameter) for treatment data
data.pred <- read.csv("../data/tree_initial_diameter.csv")
data.pred <- subset(data.pred, stem_id==1) #only use 1 stem per tree to maintain independence
newx <- data.pred$diameter_cm
newspp <- data.pred$species_name
biomass.2010 <- as.numeric(predict(model, newdata=data.frame(x=newx,species=newspp),
                                   interval="prediction", type="response"))
data.pred$initial_biomass_g <- round(exp(biomass.2010),0)
write.csv(data.pred, "../data/tree_initial_biomass.csv", row.names=FALSE)

# Pull in data to predict biomass chnage for year 1 (from dendrometer bands)
data.pred <- read.csv("../data/control_trees_growth_2011-2013.csv")
data.pred <- subset(data.pred, Stem_ID==1) #only use 1 stem per tree to maintain independence
# First estimate diameters from circumference changes for each year
data.pred$initial_diameter <- data.pred$Initial_Circumference_cm/pi
data.pred$year1_diameter <- with(data.pred, initial_diameter+(Y1_Circumference_change_cm/pi))
data.pred$year3_diameter <- with(data.pred, initial_diameter+(Y3_Circumference_change_cm/pi))

newx <- data.pred$initial_diameter
newspp <- data.pred$Species
biomass.2010 <- as.numeric(predict(model, newdata=data.frame(x=newx,species=newspp),
                                   interval="prediction", type="response"))

newx <- data.pred$year1_diameter
biomass.2011 <- as.numeric(predict(model, newdata=data.frame(x=newx,species=newspp),
                                   interval="prediction", type="response"))

newx <- data.pred$year3_diameter
newx[which(is.na(newx)==TRUE)] <- 0
biomass.2013 <- as.numeric(predict(model, newdata=data.frame(x=newx,species=newspp), 
                                   interval="prediction", type="response"))
biomass.2013[which(biomass.2013==-Inf)] <- NA

dataControlBiomass <- data.frame(Site = data.pred$Site,
                                 TreeID = data.pred$Tree_ID,
                                 Species = data.pred$Species,
                                 B2010 = round(exp(biomass.2010),0),
                                 B2011 = round(exp(biomass.2011),0),
                                 B2013 = round(exp(biomass.2013),0))

dataControlBiomass$G2011 <- with(dataControlBiomass, log(B2011)-log(B2010))
dataControlBiomass$G2013 <- with(dataControlBiomass, (log(B2013)-log(B2010))/3)

dM <- melt(dataControlBiomass, id.vars = c("Site", "TreeID", "Species"), 
           measure.vars = c("G2011", "G2013"))
dM$Year <- rep(c("2011", "2013"), each=nrow(dataControlBiomass))
dM <- subset(dM, Site != "Tiorola")
colnames(dM)[5] <- "rgr"

####
####  Write predictions to file ----------------------------
####
write.csv(dM, "../results/control_predicted_biomass.csv")

