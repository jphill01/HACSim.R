### Run HAC Simulations ###

## Input parameters ###

# Required #

# N = Number of specimens (DNA sequences)
# Hstar = Number of observed unique haplotypes
# probs = Probability frequency distribution of haplotypes
# perms = Number of permutations (replications)
# p = Proportion of unique haplotypes to recover

# Optional #

# input.seqs = Analyze inputted aligned/trimmed FASTA DNA sequence file (TRUE / FALSE)?
# subset.seqs = Take random subsample of DNA sequences (TRUE/FALSE)?
# prop.seqs = Proportion of DNA sequences to sample 
# prop.haps = Proportion of haplotypes to sample 
# subset.haps = Random subsample of haplotypes to analyze

# Run algorithm with N = 10, 50, 100, with prespecified H*, probs and p

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

N <- 50 # total number of sampled individuals
Hstar <- 10  # total number of haplotypes
probs <- c(0.45, 0.45, rep(0.1/8, 8)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency

perms <- 10000 # number of permutations
p <- 0.95 # proportion of haplotypes to recover

## Simulate hypothetical species WITHOUT migration/gene flow ##

input.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
subset.haps <- NULL # subset haplotypes?
prop.haps <- NULL # proportion of haplotypes to subsample

## Simulate hypothetical species WITH migration/gene flow ##

input.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
prop.haps <- 0.20 # proportion of haplotypes to subsample

if (!is.null(prop.haps)) { # take random subsample of haplotypes for hypothetical species
  subset.haps <- sort(sample(Hstar, size = ceiling(prop.haps * Hstar), replace = FALSE))
}

## Simulate real species WITHOUT migration/gene flow ##

input.seqs <- TRUE # analyze DNA sequence file?
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE

## Simulate real species WITH migration/gene flow ##

input.seqs <- TRUE # analyze DNA sequence file?
subset.seqs <- TRUE # subset DNA sequences?
prop.seqs <- 0.20 # proportion of DNA sequences to subsample
subset.haps <- NULL # subset haplotypes?  DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE


##########

### Run simulations ###

HAC.simrep()

##########

# Models

k <- 40
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

