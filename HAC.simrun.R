### Run HAC Simulations ###

# Run algorithm with N = 10, 50, 100, 500, 1000, with prespecified H*, probs K and p

##########

##### Set working directory #####

setwd("/Users/jarrettphillips/Desktop/HAC simulation")

##### Clear memory #####

remove(list = ls())

##### Load packages #####

library(Rcpp)
library(RcppArmadillo)
library(pegas)

### Load scripts ###

sourceCpp("accumulate.cpp")
source("HAC.sim.R")
source("HAC.simrep.R")

library(boot) # This package is for bootstrapping
library(investr) # This package performs inverse estimation
library(rootSolve) # This package employs bisection and Newton's method
library(mgcv) # This package fits GAMs and Kriging models
library(scam) # This package fits SCAMs
library(snow) # This package is for parallel processing

source("HAC.simmodels.R")
source("HAC.simplot.R")
source("HAC.simest.R")
source("HAC.simfit.R")
source("HAC.simboot.R")
source("inv.predict.R")

### Set parameters ###

# Simulate hypothetical species

N <- 100 # total number of sampled individuals
Hstar <- 10  # total number of haplotypes
probs <- c(0.45, 0.45, rep(0.10/8, 8))
# probs <- rep(1/Hstar, Hstar)
K <- 1 # number of equally-sized (sub)populations
perms <- 10000 # number of permutations
p <- 0.95 # proportion of haplotypes to recover
input.seqs <- FALSE

# Simulate real species

K <- 1 # number of equally-sized (sub)populations
perms <- 10000 # number of permutations
p <- 0.80 # proportion of haplotypes to recover
input.seqs <- TRUE

##########

### Run simulations ###

ptm <- proc.time()

HAC.simrep()

proc.time() - ptm


##########

# Models

HAC.simmodels(k = 40)

# Visualization plots

HAC.simplot(model = "GAM", k = 40)
HAC.simplot(model = "SCAM", k = 40)
HAC.simplot(model = "Krig", k = 40)

# AIC

HAC.simfit(model = "GAM", k = 40)
HAC.simfit(model = "SCAM", k = 40)
HAC.simfit(model = "Krig", k = 40)

# Parameter Estimation

HAC.simest(model = "GAM", k = 40)
HAC.simest(model = "SCAM", k = 40)
HAC.simest(model = "Krig", k = 40)

# Bootstrap simulation

HAC.simboot(model = "GAM", k = 40)
HAC.simboot(model = "SCAM", k = 40)
HAC.simboot(model = "Krig", k = 40)

