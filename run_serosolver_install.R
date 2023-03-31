remotes::install_github("seroanalytics/serosolver", ref = "published", 
                        dependencies = TRUE,force = TRUE)
#devtools::install_github("seroanalytics/serosolver") # this one definately does not work.
library(serosolver)

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
