### Run HAC Simulations ###

## Input parameters ###

# Required #

# N = Number of specimens (DNA sequences)
# Hstar = Number of observed unique haplotypes
# probs = Probability frequency distribution of haplotypes

# Optional #

# p = Proportion of unique haplotypes to recover
# perms = Number of permutations (replications) 
# input.seqs = Analyze inputted aligned/trimmed FASTA DNA sequence file (TRUE / FALSE)?
# sim.seqs = Simulate DNA sequences (TRUE / FALSE)?
# num.seqs = Number of DNA sequences to simulate
# length.seqs = Basepair length of DNA sequences to simulate
# mu.rate = Substitution rate of simulated DNA sequences
# subst.model = Nucleotide substition model
# transi.rate = Substitution rate of transitions of simulated DNA sequences under K2P model
# transv.rate = Substitution rate of transversions of simulated DNA sequences under K2P model
# subset.seqs = Subset of DNA sequences to sample
# prop.seqs = Proportion of DNA sequences to sample 
# prop.haps = Proportion of haplotypes to sample 
# subset.haps = Subset of haplotypes to sample
# num.pts = Number of terminal data points used to calculate curve slope 
# prop.pts = Proportion of terminal data points used to calculate curve slope


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
#source("HAC.simbestaic.R")
#source("HAC.simplot.R")
#source("HAC.simest.R")
#source("HAC.simfit.R")
#source("HAC.simmodels.R")
#source("HAC.simboot.bisect.R")
#source("HAC.simboot.newton.R")
#source("HAC.simpost.R")
#source("inv.predict.R")



##### SCENARIOS #####

# N = 10, 50, 100
# Hstar = 10, 20, 25

## Example probs for Hstar = 10 - you can come up with others ##

# probs = rep(1/Hstar, Hstar) # equal haplotype frequency
# probs = c(0.90, rep(0.10/9, 9)) # 1 dominant haplotype
# probs = c(0.45, 0.45, rep(0.10/8, 8)) # 2 dominant haplotypes
# probs = c(0.30, 0.30, 0.30, rep(0.10/7, 7)) # 3 dominant haplotypes

## Haplotype Recovery ## 

# p = 0.90
# p = 0.95 - default
# p = 0.99
# p = 1 - long runtime


### SET PARAMETERS ###

## Run for all simulations ##

perms <- 10000 # number of permutations
p <- 0.95 # proportion of haplotypes to recover
num.pts <- 10 # number of terminal data points for curve slope calculation
prop.pts <- NULL # proportion of terminal data points for curve slope calculation


## Set parameters for hypothetical species ##

N <- 100 # total number of sampled individuals
Hstar <- 15 # total number of haplotypes
# probs <- c(53/76, 13/76, 9/76, 1/76) # must sum to 1
probs <- c(214/234, rep(3/234, 2), rep(2/234, 2), rep(1/234, 10))
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency


## Simulate hypothetical species WITHOUT migration/gene flow ##

input.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
sim.seqs <- FALSE # simulate DNA sequrnces? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE

## Simulate hypothetical species WITH migration/gene flow ##

input.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
sim.seqs <- FALSE # simulate DNA sequrnces? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
prop.haps <- 0.50 # proportion of haplotypes to subsample

# subset.haps cannot have a length of 1

if (!is.null(prop.haps)) { # take random subsample of haplotypes for hypothetical species
  subset.haps <- sort(sample(Hstar, size = ceiling(prop.haps * Hstar), replace = FALSE))
}


## Simulate real species WITHOUT migration/gene flow ##

input.seqs <- TRUE # analyze DNA sequence file? DO NOT CHANGE
sim.seqs <- FALSE # simulate DNA sequrnces? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE

## Simulate real species WITH migration/gene flow ##

input.seqs <- TRUE # analyze DNA sequence file? DO NOT CHANGE
sim.seqs <- FALSE # simulate DNA sequrnces? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes?  DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- TRUE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- 0.10 # proportion of DNA sequences to subsample


## Simulate DNA sequences ## 

# JC69: mu.transi = mu.transv, A = C = G = T = 0.25
# K80: mu.transi != mu.transv, A = C = G = T = 0.25

# F81: mu.transi = mu.transv, A != C != G != T != 0.25
# HKY85: mu.transi != mu.transv, A != C != G != T != 0.25


# JC69

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE 
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
num.seqs <- 100 # number of DNA sequences to simulate
length.seqs <- 658 # length of DNA sequences to simulate
subst.model <- "JC69" # nucleotide substitution model DO NOT CHANGE
nucl.freq <- NULL # nucleotide frequencies DO NOT CHANGE
mu.rate <- 1e-4 # mutation rate
transi.rate <- NULL # transition rate DO NOT CHANGE
transv.rate <- NULL  # transversion rate DO NOT CHANGE


# K80

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
num.seqs <- 100 # number of DNA sequences to simulate
length.seqs <- 658 # length of DNA sequences to simulate
subst.model <- "K80" # nucleotide substitution model DO NOT CHANGE
nucl.freq <- NULL # nucleotide frequencies DO NOT CHANGE
mu.rate <- NULL # mutation rate DO NOT CHANGE
transi.rate <- 1e-4 # transition rate
transv.rate <- transi.rate / 2  # transversion rate


# F81

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
num.seqs <- 100 # number of DNA sequences to simulate
length.seqs <- 658 # length of DNA sequences to simulate
subst.model <- "F81" # nucleotide substitution model DO NOT CHANGE
nucl.freq <- c(0.70, 0.10, 0.10, 0.10)
mu.rate <- 1e-4 # mutation rate
transi.rate <- NULL # transition rate DO NOT CHANGE
transv.rate <- NULL  # transversion rate DO NOT CHANGE


# HKY85

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
num.seqs <- 100 # number of DNA sequences to simulate
length.seqs <- 658 # length of DNA sequences to simulate
subst.model <- "HKY85" # nucleotide substitution model DO NOT CHANGE
nucl.freq <- c(0.70, 0.10, 0.10, 0.10)
mu.rate <- NULL # mutation rate DO NOT CHANGE
transi.rate <- 1e-4 # transition rate
transv.rate <- transi.rate / 2  # transversion rate



### Run simulations ###

HAC.simrep(filename = "output")


##########

# Models

k <- 10 # default k - likely will have to double each time until k is big enough
HAC.simmodels(k = k)

# Visualization plots - check to 

HAC.simfit(model = "GAM", k = k) # check to ensure k is big enough - dashed line should be monotone increasing and go through all data points
HAC.simfit(model = "SCAM", k = k)
HAC.simfit(model = "Krig", k = k)

HAC.simplot(model = "GAM", k = k) # histograms, QQplots - check to ensure k is big enough
HAC.simplot(model = "SCAM", k = k)
HAC.simplot(model = "Krig", k = k)

# AIC for all 10 models

HAC.simaic(model = "GAM", k = k)
HAC.simaic(model = "SCAM", k = k)
HAC.simaic(model = "Krig", k = k)

# Best Model + AIC

HAC.simbestaic(model = "GAM")
HAC.simbestaic(model = "SCAM")
HAC.simbestaic(model = "Krig")
HAC.simbestaic(model = "All")

# Bootstrap simulation - CAN BE SLOW 

HAC.simboot.bisect(model = "GAM", k = k)
HAC.simboot.bisect(model = "SCAM", k = k)
HAC.simboot.bisect(model = "Krig", k = k)

HAC.simboot.newton(model = "GAM", k = k)
HAC.simboot.newton(model = "SCAM", k = k)
HAC.simboot.newton(model = "Krig", k = k)

