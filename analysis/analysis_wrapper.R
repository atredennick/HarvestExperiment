##  Wrapper script to call individual R analysis scripts
##    to completely reproduce all results and figures from
##    Tredennick et al. 2015, African Journal of Ecology.

##  Run this script, and all results and figures will be found
##    in the results/ folder, which is made on the fly.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4.7.2015
##                4.20.2015 -- Adds check for recent Rcpp package.
##                          -- Adds Rcpp to packages list.

####
#### Install and load necessary packages -----------------------------------------
####
packages <- c("ggplot2", "plyr", "reshape2", "rdryad", "Rcpp",
              "ggthemes", "lme4", "gridExtra", "AICcmodavg")
# Uncomment next line if needed...
# install.packages(packages)

# Check for packages
for(pack in 1:length(packages)){
  if(packages[pack] %in% rownames(installed.packages())==FALSE)
  {stop(paste("You need to install the", packages[pack]), " package from CRAN.")}
}

# Auto update Rcpp package if not 0.11.5
# This will update even if package loaded is beyond 0.11.5
# This is just a failsafe to make sure Rcpp will work
if(packageVersion("Rcpp") != "0.11.5")
  install.packages("Rcpp")

####
####  Create new results/ directory and run analysis -----------------------------
####
dir.create("../results/", showWarnings = FALSE)

# Run all analyses, in this order...
source("get_allometric_data.R")
source("allometry_linear_model.R")
source("control_fh_ANOVA.R")
source("treatment_ANOVA.R")
source("height_analysis.R")
source("wet_dry_ratio_test.R")

print("All analysis for Tredennick et al. 2015 complete. Check your results folder.")
print("If errors occur or you have questions, email Andrew Tredennick at atredenn@gmail.com.")


