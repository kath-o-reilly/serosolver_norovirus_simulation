## Read in chains for trace plot
# requires
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(readr)
library(viridis)
library(ggpubr)
library(serosolver)

## MCMC settings, not super important but can be tweaked
mcmc_pars <- c("save_block"=100,"thin"=100,"thin_hist"=500,"iterations"=100000,
               "adaptive_period"=50000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)


tmp1 <- "/Users/lsh1603970/GitHub/serosolver_norovirus_simulation/chains/sim_noro" # works fine
chains <- load_mcmc_chains(tmp1,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)

rm(chains)

tmp2 <- "/Users/lsh1603970/GitHub/serosolver_norovirus_simulation/chains/sim_noro_allsampled" # works fine
chains <- load_mcmc_chains(tmp2,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
