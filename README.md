No effects of fire, large herbivores, and their interaction on regrowth of harvested trees in two West African savannas
=================

Data and code to reproduce analysis of SSDE harvest experiment as reported in Tredennick et al. (in review).

The abstract of our paper:
>Theory and empirical evidence for the impacts of fire and herbivory in savannas is well established – they are top-down disturbances that maintain savannas in disequilibrium states away from potential tree cover. In African savannas the demand for fuelwood is extremely high, so tree harvest likely also has an impact, both directly and indirectly, on tree cover, density, and biomass. Many savanna trees resprout vigorously from the base after harvest. However, harvested trees regenerate as saplings susceptible to fire and browsing, so harvest may have important demographic consequences. Here, we report the effects of tree harvest, and its interaction with fire and herbivory, on savanna dynamics by analyzing woody regrowth following a harvest in arid Sahelian and mesic Guinean savannas in Mali, West Africa. Tree harvest resulted in an overall reduction in wood production per tree compared to growth in non-harvested trees. Regrowth, either biomass or height, did not differ among fire and herbivory treatments. Our results suggest that the resprouting abilities that savanna trees have evolved to cope with frequent fire are essential for surviving tree harvest and subsequent disturbance. In these savannas, regrowth is rapid enough in the first growing season to escape the impact of dry season fires.

Requirements
------------------------
For all of our results, we used R (http://www.r-project.org/; we used R 3.1.0, but our code should be compatible across builds) and the following R packges:

* ``gridExtra``
* ``plyr``
* ``reshape2``
* ``ggplot2``
* ``lme4``
* ``rdryad``
* ``ggthemes``
* ``AICcmodavg``

To install these packages, run the code chunk below in your R console:

```coffee
install.packages("ggplot2", dependencies = TRUE)
install.packages("plyr")
install.packages("reshape2")
install.packages("gridExtra")
install.packages("lme4")
install.packages("rdryad")
install.packages("AICcmodavg")
install.packages("ggthemes")
```
