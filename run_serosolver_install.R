
# option 1 (default)

remotes::install_github("seroanalytics/serosolver", ref = "published", 
                        dependencies = TRUE,force = TRUE)
#devtools::install_github("seroanalytics/serosolver") # this one definately does not work.
library(serosolver)

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# option 2 (2 sera option)

# remove.packages(serosolver) # *** this is only for things that have been loaded?

# *** seems to work? -- use the packages tab and remove from here --
# Removing package from ‘/Library/Frameworks/R.framework/Versions/4.2/Resources/library’
#(as ‘lib’ is unspecified)

devtools::install_github("seroanalytics/serosolver",ref="multiple_obs_types") # installs the package directly from github
