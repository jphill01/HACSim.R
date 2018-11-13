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
# subst.model = Nucleotide substition model (JC69 / K80 / F81 / HKY85)
# subset.seqs = Subset of DNA sequences to sample
# prop.seqs = Proportion of DNA sequences to sample 
# prop.haps = Proportion of haplotypes to sample 
# subset.haps = Subset of haplotypes to sample
# num.pts = Number of points used to calculate curve slope 
# prop.pts = Proportion of points used to calculate curve slope

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
# library(pegas)
# library(rootSolve) # This package employs bisection and Newton's method
# library(mgcv) # This package fits GAMs and Kriging models
# library(scam) # This package fits SCAMs

# source("HAC.simmodels.R")
# source("HAC.simaic.R")
# source("HAC.simbestaic.R")
# source("HAC.simplot.R")
# source("HAC.simfit.R")
# source("inv.predict.R")

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

## Set parameters for hypothetical species ##

N <- 100 # total number of sampled individuals
Hstar <- 10 # total number of haplotypes
probs <- c(rep(0.30, 3), rep(0.10/7, 7)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency


## Simulate hypothetical species WITHOUT migration/gene flow ##

input.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
sim.seqs <- FALSE # simulate DNA sequences? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE


## Simulate hypothetical species WITH migration/gene flow ##

input.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
sim.seqs <- FALSE # simulate DNA sequrnces? DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
prop.haps <- 0.15 # proportion of haplotypes to subsample

# Subsample haplotypes ##

subset.haps <- sort(sample(Hstar, size = ceiling(prop.haps * Hstar), replace = FALSE))


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

# raw: proportion of nucleotide sites that differ among specimens

# JC69: corrected p-distance; mu.transi = mu.transv, A = C = G = T = 0.25
# K80: mu.transi != mu.transv, A = C = G = T = 0.25

# F81: mu.transi = mu.transv, A != C != G != T != 0.25


# p-distance

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
subst.model <- "raw" # nucleotide substitution model DO NOT CHANGE

# JC69

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE 
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
subst.model <- "JC69" # nucleotide substitution model DO NOT CHANGE

# K80

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
subst.model <- "K80" # nucleotide substitution model DO NOT CHANGE

# F81

input.seqs <- FALSE # analyze DNA sequence file? DO NOT CHANGE
subset.haps <- NULL # subset haplotypes? DO NOT CHANGE
prop.haps <- NULL # proportion of haplotypes to subsample DO NOT CHANGE
subset.seqs <- FALSE # subset DNA sequences? DO NOT CHANGE
prop.seqs <- NULL # proportion of DNA sequences to subsample DO NOT CHANGE
sim.seqs <- TRUE # simulate DNA sequences? DO NOT CHANGE
subst.model <- "F81" # nucleotide substitution model DO NOT CHANGE


### Run simulations ###

HAC.simrep(filename = "output")


##########

# Models

k <- 10 # default k - likely will have to double each time until k is big enough
HAC.simmodels(k = k)

# Visualization plots 

HAC.simfit(model = "GAM", k = k) # check to ensure k is big enough - dashed line should be monotone increasing and go through all data points
HAC.simfit(model = "SCAM", k = k)
HAC.simfit(model = "Krig", k = k)

HAC.simplot(model = "GAM", k = k) # histograms, QQplots - check to ensure k is big enough
HAC.simplot(model = "SCAM", k = k)
HAC.simplot(model = "Krig", k = k)


# Model AICs

HAC.simaic(model = "GAM")
HAC.simaic(model = "SCAM")
HAC.simaic(model = "Krig")
HAC.simaic(model = "Best")
HAC.simaic(model = "All")


## Large sample (Asymptotically Normal ) CIs ##

HAC.simest(model = "GAM", k = k, method = "Bisect")
HAC.simest(model = "SCAM", k = k, method = "Bisect")
HAC.simest(model = "Krig", k = k, method = "Bisect")

HAC.simest(model = "GAM", k = k, method = "Newton")
HAC.simest(model = "SCAM", k = k, method = "Newton")
HAC.simest(model = "Krig", k = k, method = "Newton")

HAC.simest(model = "GAM", k = k, method = "Both")
HAC.simest(model = "SCAM", k = k, method = "Both")
HAC.simest(model = "Krig", k = k, method = "Both")


## Bootstrap simulation - 1000 reps - CAN BE SLOW ##

# Individual models

HAC.simboot(model = "GAM", k = k, bootType = "Bisect", boot.reps = 1000)
HAC.simboot(model = "SCAM", k = k, bootType = "Bisect", boot.reps = 1000)
HAC.simboot(model = "Krig", k = k, bootType = "Bisect", boot.reps = 1000)

HAC.simboot(model = "GAM", k = k, bootType = "Newton", boot.reps = 1000)
HAC.simboot(model = "SCAM", k = k, bootType = "Newton", boot.reps = 1000)
HAC.simboot(model = "Krig", k = k, bootType = "Newton", boot.reps = 1000)

# Best model

HAC.simboot(model = "Best", k = k, bootType = "Bisect", boot.reps = 1000)
HAC.simboot(model = "Best", k = k, bootType = "Newton", boot.reps = 1000)

# All models

HAC.simboot(model = "All", k = k, bootType = "Bisect", boot.reps = 1000)
HAC.simboot(model = "All", k = k, bootType = "Newton", boot.reps = 1000)
