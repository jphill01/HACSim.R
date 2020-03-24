##### Set working directory #####

setwd("/Users/jarrettphillips/Desktop")

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

N <- 76 # total number of sampled individuals
Hstar <- 4 # total number of haplotypes
# probs <- c(rep(0.45, 2), rep(0.10/4, 4)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency
probs <- c(53, 13, 9, 1) / N

### Run HAC Simulations ###

## Simulate hypothetical species ##
## If not set, prop defaults to 0.10
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 1000, p = 0.95, conf.level = 0.95, num.iters = NULL, filename = "output")

## Simulate hypothetical species - subsampling ##
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 1000, p = 0.95, subsample = TRUE, prop = 0.75, conf.level = 0.95, num.iters = NULL, filename = "output")

## Simulate hypothetical species and all paramaters changed - subsampling ##
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.95, filename = "output")

## Simulate real species ##
HACSObj <- HACReal(p = 0.95, perms = 10000, conf.level = 0.95, num.iters = NULL, filename = "output")

## Simulate real species - subsampling ##
HACSObj <- HACReal(p = 0.95, subsample = TRUE, prop = 0.25, conf.level = 0.95, num.iters = NULL, filename = "output")

## Simulate real species and all parameters changed - subsampling ##
HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.99, filename = "output")

## Simulate DNA sequences

sim.seqs(num.seqs = 100, num.haps = 8, length.seqs = 658, nucl.freqs = rep(0.25, 4), count.haps = c(50, 20, rep(5, 6)), subst.model = "JC69", mu.rate = 1e-4)

HACSObj <- HACReal(p = 0.95, perms = 10000, conf.level = 0.95, num.iters = NULL, filename = "output")

set.seed(0673227) # set random seed for reproducibility, if desired
HAC.simrep(HACSObj) # simulate
