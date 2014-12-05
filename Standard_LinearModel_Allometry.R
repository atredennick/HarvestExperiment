rm(list=ls()) 
library(ggplot2)
library(nlme)
library(plyr)
library(reshape2)

###NOTE, THIS MODEL PREDICTS THE LOG OF BIOMASS IN GRAMS###

####### WHOLE TREE-LEVEL #######
##Bring in data
data2=read.csv("Allom_Data_MS.csv")
names(data2) = tolower(names(data2))

data2 = data2[data2$aggregated_leaf_wt>-9000,]

Y2 = data2$aggregated_wood_wt
x1.2 = data2$diameter
x2.2 = data2$length
species = data2$species


## TEST MODEL WITH DIAMETER ONLY ##
###### MIXED-EFFECTS MODEL -- MAIN STEMS ########
## Species as random effect 
species_null = (numeric(length(species)))

##With species random effect
mix1 = lme(fixed=log(Y2)~log(x1.2), random = ~ 1 | species, method="ML")
summary(mix1)

##Update to remove species effect
mix2 = update(mix1, random = ~ 1 | species_null)
summary(mix2)

##Compare models
anova(mix1, mix2)
AICc = function(AIC, n, k){
	AICc = AIC + 2*k * (n/(n-k-1))
	return(AICc)
	}

AICc1=AICc(AIC(mix1), n=38, k=3)
AICc2=AICc(AIC(mix2), n=38, k=2)

## Final Model ##
x=x1.2
model <- lm(log(Y2)~log(x))
summary(model)

#quickplot
qplot(log(x),log(Y2), size=2)+
  geom_line(aes(x=log(x), y=predict(model)), size=1.5)+
  xlab("log(Diameter)")+
  ylab("log(Wood Biomass)")+
  guides(size=FALSE)

##Pull in data to predict biomass chnage for year 1 -- from dendrometer bands
data.pred2 <- read.csv("Control Trees Growth_2011-2013.csv")

newx=data.pred2$Diameter..year.1.
biomass.2010 <- predict(model, newdata=data.frame(x=newx), interval="prediction", type="response")

newx=data.pred2$D2_Y2..cm.
biomass.2011 <- predict(model, newdata=data.frame(x=newx), interval="prediction", type="response")

newx=data.pred2$D3_Y3
biomass.2013 <- predict(model, newdata=data.frame(x=newx), interval="prediction", type="response")

dataControlBiomass <- data.frame(Site = data.pred2$Site,
                                 TreeID = data.pred2$Tree_ID,
                                 Species = data.pred2$Species,
                                 B2010 = exp(biomass.2010[,1]),
                                 B2011 = exp(biomass.2011[,1]),
                                 B2013 = exp(biomass.2013[,1]))

dataControlBiomass$G2011 <- with(dataControlBiomass, log(B2011)-log(B2010))
dataControlBiomass$G2013 <- with(dataControlBiomass, (log(B2013)-log(B2010))/3)

dM <- melt(dataControlBiomass, id.vars = c("Site", "TreeID", "Species"), 
           measure.vars = c("G2011", "G2013"))
dM$Year <- rep(c("2011", "2013"), each=nrow(dataControlBiomass))
dM <- subset(dM, Site != "Tiorola")
colnames(dM)[5] <- "rgr"

write.csv(dM, "controlPredBiomass.csv")

ggplot(data=dM, aes(x=Site, y=rgr, fill=Year))+
  geom_boxplot(outlier.shape = 1)+
  ylab(expression(paste("Relative Growth Rate (g ", g^-1, " ", y^-1,")")))+
  scale_fill_manual(labels=c("2011", "2013"), values=c("grey40", "grey80"),
                    name="Year")


