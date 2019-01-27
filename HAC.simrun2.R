##### Set working directory #####

setwd("/Users/jarrettphillips/Desktop/HAC simulation")

##### Clear memory #####

remove(list = ls())

##### Load packages #####

library(HACSim)
?HACSim


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

N <- 30 # total number of sampled individuals
Hstar <- 10 # total number of haplotypes
probs <- c(rep(0.45, 2), rep(0.10/8, 8)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency

### Run HAC Simulations ###

## Simulate hypothetical species ##
## If not set, prop defaults to 0.10
HACSObj <- HACHypothetical(N, Hstar, probs, filename = "output")

## Simulate hypothetical species - subsampling ##
HACSObj <- HACHypothetical(N, Hstar, probs, subsample = TRUE, prop = 0.25, filename = "output")

## Simulate hypothetical species and all paramaters changed - subsampling ##
HACSObj <- HACHypothetical(N, Hstar, probs, perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, filename = "output")

## Simulate real species ##
HACSObj <- HAC(type = "real", filename = "output")

## Simulate real species - subsampling ##
HACSObj <- HAC(type = "real", subsample = TRUE, prop = 0.15, filename = "output")

## Simulate real species and all paramaters changed - subsampling ##
HACSObj <- HAC(type = "real", perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, filename = "output")

# Simulate DNA sequences under JC69
HACSObj <- HAC(type = "evolution", subst.model = "JC69", num.seqs = 100, length.seqs = 500, mu.rate = 0.0010)

HAC.simrep(HACSObj)
