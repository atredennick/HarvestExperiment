##  This script downloads the original allometric data from
##    from Dryad and pulls the main stem records we want for
##    estimating the whole-tree allometric model for relating
##    tree diameter to tree biomass.

##  This script requires and internet connection.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.8.2015

library('rdryad')

allom_data <- dryad_getfile("http://datadryad.org/bitstream/handle/10255/dryad.46421/Mali_Savanna_AllometricData.csv?sequence=1")
allom_data <- subset(allom_data, Site!="Tiorola")
main_stem_records <- grep("MS", allom_data$Branch_ID)
mains_data <- allom_data[main_stem_records,]
write.csv(mains_data, "../data/allometry_data_main_stems.csv")
