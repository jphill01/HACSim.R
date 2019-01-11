##### Set working directory #####

setwd("/Users/jarrettphillips/Desktop/HAC simulation")

##### Clear memory #####

remove(list = ls())

##### Load packages #####

library(HACSim)
?HACSim

### Load scripts ###

# library(Rcpp)
# library(RcppArmadillo)

# sourceCpp("accumulate.cpp")
# source("HAC.sim.R")
# source("HAC.simrep.R")

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

## Set parameters for hypothetical species ##

N <- 100 # total number of sampled individuals
Hstar <- 10 # total number of haplotypes
probs <- c(rep(0.45, 2), rep(0.10/8, 8)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency

### Run HAC Simulations ###

## Simulate hypothetical species WITHOUT migration/gene flow ##
## If not set, prop defaults to 0.10
HACSObj <- HACHypothetical(N, Hstar, probs)

## Simulate hypothetical species WITH migration/gene flow ##
HACSObj <- HACHypothetical(N, Hstar, probs, subsample = TRUE, prop = 0.15)

## Simulate hypothetical species WITH migration/gene flow and all paramaters changed ##
HACSObj <- HACHypothetical(N, Hstar, probs, perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, filename = "output")

## Simulate real species WITHOUT migration/gene flow ##
HACSObj <- HAC(type = "real")

## Simulate real species WITH migration/gene flow ##
HACSObj <- HAC(type = "real", subsample = TRUE, prop = 0.15)

## Simulate real species WITH migration/gene flow and all paramaters changed, these can be changes for evolution as well##
HACSObj <- HAC(type = "real", perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, filename = "output")

## Simulate DNA sequences ## 

# p-distance
HACSObj <- HAC(type = "evolution", subst.model = "raw")

# JC69
HACSObj <- HAC(type = "evolution",subst.model = "JC69")

# K80
HACSObj <- HAC(type = "evolution", subst.model = "K80")

# F81 
HACSObj <- HAC(type = "evolution", subst.model = "F81")

# JC69 and all parameters changed ##
HACSObj <- HAC(type = "evolution", perms = 10000, p = 0.90, subsample = TRUE, subst.model = "JC69", filename = "output")


HAC.simrep(HACSObj)


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
