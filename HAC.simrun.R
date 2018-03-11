### Run HAC Simulations ###

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
# sourceCpp("migrate.cpp")
source("HAC.sim.R")
source("HAC.simrep.R")

library(boot) # This package is for bootstrapping
library(investr) # This package performs inverse estimation
library(rootSolve) # This package employs the bisection method
library(mgcv) # This package fits GAMs and Kriging models
library(scam) # This package fits SCAMs
library(doParallel)

source("HAC.simplot.R")
source("HAC.simest.R")
source("HAC.simfit.R")
source("HAC.simboot.R")
source("inv.predict.R")

### Set parameters ###

# Simulate hypothetical species

N <- 10 # total number of sampled individuals
Hstar <- 5  # total number of haplotypes
#probs <- c(0.30, 0.30, 0.30, rep(0.10/2, 2))
probs <- rep(1/Hstar, Hstar)
K <- 1 # number of equally-sized (sub)populations
perms <- 10000 # number of permutations
p <- 0.6 # proportion of haplotypes to recover
input.seqs <- FALSE

# Simulate real species

K <- 1 # number of equally-sized (sub)populations
perms <- 10000 # number of permutations
p <- 1 # proportion of haplotypes to recover
input.seqs <- TRUE

##########

### Run simulations ###

HAC.sim(N = N, Hstar = Hstar, probs = probs, K = K, perms = perms, p = p, input.seqs = input.seqs)
HAC.simrep()

##########

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