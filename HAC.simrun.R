### Run HAC Simulations ###

# Run algorithm with N = 10, 50, 100, with prespecified H*, probs, K, m and p

##########

##### Set working directory #####

setwd("/Users/jarrettphillips/Desktop/HAC simulation")

##### Clear memory #####

remove(list = ls())

##### Load packages #####

library(HACSim)
?HACSim

### Load scripts ###

#library(Rcpp)
#library(RcppArmadillo)

#sourceCpp("accumulate.cpp")
#source("HAC.sim.R")
#source("HAC.simrep.R")

# library(boot) # This package is for bootstrapping
# library(investr) # This package performs inverse estimation
# library(rootSolve) # This package employs bisection and Newton's method
# library(mgcv) # This package fits GAMs and Kriging models
# library(scam) # This package fits SCAMs

#source("HAC.simmodels.R")
#source("HAC.simaic.R")
#source("HAC.simplot.R")
#source("HAC.simest.R")
#source("HAC.simfit.R")
#source("HAC.simmodels.R")
#source("HAC.simboot.bisect.R")
#source("HAC.simboot.newton.R")
#source("HAC.simpost.R")
#source("inv.predict.R")

### Set parameters ###

# Simulate hypothetical species

N <- 10 # total number of sampled individuals
Hstar <- 5  # total number of haplotypes
probs <- c(0.45, 0.45, rep(0.10/3, 3)) # must sum to 1
# probs <- rep(1/Hstar, Hstar)
perms <- 100000 # number of permutations
p <- 0.95 # proportion of haplotypes to recover
input.seqs <- FALSE

# Simulate real species

perms <- 10000 # number of permutations
p <- 0.80 # proportion of haplotypes to recover
input.seqs <- TRUE


##########

### Run simulations ###
HAC.simrep()


##########

# Models

k <- 10
HAC.simmodels(k = k)

# Visualization plots

HAC.simfit(model = "GAM", k = k)
HAC.simfit(model = "SCAM", k = k)
HAC.simfit(model = "Krig", k = k)

HAC.simplot(model = "GAM", k = k)
HAC.simplot(model = "SCAM", k = k)
HAC.simplot(model = "Krig", k = k)

# AIC

HAC.simaic(model = "GAM", k = k)
HAC.simaic(model = "SCAM", k = k)
HAC.simaic(model = "Krig", k = k)

# Parameter Estimation - CI is outputted 

HAC.simest(model = "GAM", k = k)
HAC.simest(model = "SCAM", k = k)
HAC.simest(model = "Krig", k = k)

# Posterior Simulation (MVN distribution)

HAC.simpost(model = "GAM", k = k)
HAC.simpost(model = "SCAM", k = k)
HAC.simpost(model = "Krig", k = k)

# Bootstrap simulation - VERY SLOW 

HAC.simboot.bisect(model = "GAM", k = k)
HAC.simboot.bisect(model = "SCAM", k = k)
HAC.simboot.bisect(model = "Krig", k = k)

HAC.simboot.newton(model = "GAM", k = k)
HAC.simboot.newton(model = "SCAM", k = k)
HAC.simboot.newton(model = "Krig", k = k)

