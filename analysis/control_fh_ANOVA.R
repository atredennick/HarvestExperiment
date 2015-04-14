##  This script compares regrowth of protected harvested (fh) trees
##    and protected non-harvest (control) trees. We use simple ANOVA.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.6.2015

##  This code is provided 'as is.' Any attempt to use this
##    code for different data is done so at the users own risk.

# Clear everything
rm(list=ls()) 

####
####  Load libraries ---------------------------------
####
#install.packages(c('ggplot2','plyr'))
library('ggplot2')
library('plyr')


####
####  Read in data and format -----------------------
####
# Get fh treatment data and convert to dry weights
wet_dry_data <- read.csv("../data/wet_dry_subsamples.csv")
wet_dry_data$ratio <- with(wet_dry_data, dry_wood_weight_g / wet_wood_weight_g)
mean_wet_dry_species <- ddply(wet_dry_data, .(species_name, year), summarise,
                              mean_wet_dry = round(mean(ratio),2))
#set Lakamane 2013 wet:dry to 0.5 (see supplementary text)
mean_wet_dry_species <- rbind(mean_wet_dry_species,
                              c("Combretum glutinosum", 2013, 0.5))
treatment_data <- read.csv("../data/tree_regrowth.csv")
treatment_data <- merge(treatment_data, mean_wet_dry_species, 
                        by.x=c("species_name", "harvest_year"),
                        by.y=c("species_name", "year"))
treatment_data$wood_dry_weight_g <- with(treatment_data,
                                         (wood_wet_weight_kg*as.numeric(mean_wet_dry))*1000)
fh_data <- subset(treatment_data, treatment_code==4)

# Merge in initial biomass data
initial_biomass <- read.csv("../data/tree_initial_biomass.csv")
initial_biomass <- subset(initial_biomass, treatment_code==4)
initial_biomass <- initial_biomass[,c("site", "tree_id", "initial_biomass_g")]
fh_data <- merge(fh_data, initial_biomass, by=c("site", "tree_id"))
fh_data$rrg <- with(fh_data, 
                    (log(wood_dry_weight_g + initial_biomass_g) - log(initial_biomass_g)) / (harvest_year-2010))
fh_data$treatment <- "H"
fh_data <- fh_data[,c("site", "species_name", "harvest_year", "rrg", "treatment")]

# Now get control tree data
control_data <- read.csv("../results/control_predicted_biomass.csv")
control_data <- control_data[,c("Site", "Species", "Year", "rgr")]
control_data$Treatment <- "C"
set_to_na <- which(control_data$rgr==0)
control_data$rgr[set_to_na] <- NA

# Combine the two datasets
colnames(fh_data) <- colnames(control_data)
allD <- rbind(control_data, fh_data)
colnames(allD)[which(colnames(allD)=="rgr")] <- "rrg"

# Transform the regrowth
allD$rrgNorm <- allD$rrg
allD$rrg <- log(allD$rrg) # log transform to meet normality
                          # testing of normality assumptions not shown

# Set any -Inf to NA

####
####  Run ANOVAs by site ----------------------------
####
mod1 <- aov(rrg~as.factor(Treatment)*as.factor(Year), data=subset(allD, Site=="Tiendega"))
anova(mod1)
T1 <- TukeyHSD(mod1)
comparisons <- row.names(T1[[3]])
row.names(T1[[3]]) <- NULL 

mod2 <- aov(rrg~as.factor(Treatment)*as.factor(Year), data=subset(allD, Site=="Lakamane"))
anova(mod2)
T2 <- TukeyHSD(mod2)
row.names(T2[[3]]) <- NULL 

tukeys <- as.data.frame(rbind(T1[[3]], T2[[3]]))
tukeys$site <- rep(c("Tiendega", "Lakamane"), each=nrow(T1[[3]]))
tukeys$comparison <- rep(comparisons, times=2)
write.csv(tukeys, "../results/control_fh_tukeys.csv")


####
#### Calculate effect sizes -------------------------
####
# Effect size = log difference
avgC <- ddply(allD, .(Year, Site, Treatment), summarise,
                  avgRRGR = mean(rrgNorm, na.rm = TRUE))
lak1 <- subset(avgC, Year==2011 & Site=="Lakamane")[1,4]/subset(avgC, Year==2011 & Site=="Lakamane")[2,4]
lak2 <- subset(avgC, Year==2013 & Site=="Lakamane")[1,4]/subset(avgC, Year==2013 & Site=="Lakamane")[2,4]
tien1 <- subset(avgC, Year==2011 & Site=="Tiendega")[1,4]/subset(avgC, Year==2011 & Site=="Tiendega")[2,4]
tien2 <- subset(avgC, Year==2013 & Site=="Tiendega")[1,4]/subset(avgC, Year==2013 & Site=="Tiendega")[2,4]

# Calculate average log difference across years
mean(c(lak1, lak2))
mean(c(tien1, tien2))
out_d <- data.frame(site = c("Tiendega", "Lakamane"),
                    logdiff = c(mean(c(tien1, tien2)), mean(c(lak1, lak2))))
write.csv(out_d, "../results/control_fh_logdifferences.csv")

####
####  Make Figure 2 ------------------------------------------------
####
pdf("../results/Figure2.pdf", width=5, height=6)
print(
ggplot(allD, aes(x=as.factor(Year), y=rrgNorm, fill=Treatment))+
  geom_boxplot(outlier.shape = 1, na.rm=TRUE)+
  facet_grid(Site~.)+
  scale_fill_manual(labels=c("Control", "Harvested (fh)"), values=c("grey40", "grey80"),
                    name="")+
  ylab(expression(paste("Relative (Re)Growth Rate (g ", g^-1, " ", y^-1,")")))+
  xlab("Measurement Year")+
  theme_bw()
)
dev.off()



