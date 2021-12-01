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

N <- 76 # total number of sampled individuals
Hstar <- 4 # total number of haplotypes
# probs <- c(rep(0.45, 2), rep(0.10/4, 4)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency
probs <- c(53, 13, 9, 1) / N

### Run HAC Simulations ###

## Simulate hypothetical species ##
## If not set, prop defaults to 0.10
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 10000, p = 0.95, conf.level = 0.95, num.iters = NULL, filename = "output")

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

## Set parameters to simulate DNA sequences ##

num.seqs <- 100 # number of DNA sequences
num.haps <- 15 # number of haplotypes
length.seqs <- 658 # length of DNA sequences
count.haps <- c(60, rep(10, 2), rep(5, 2), rep(1, 10)) # haplotype frequency distribution
nucl.freqs <- c(0.4, 0.3, 0.2, 0.1) # nucleotide frequency distribution
subst.model <- "F81" # desired nucleotide substitution model
mu.rate <- 1e-3 # mutation rate
transi.rate <- NULL # transition rate
transv.rate <- NULL # transversion rate

sim.seqs(num.seqs = num.seqs, num.haps = num.haps, length.seqs = length.seqs, nucl.freqs = nucl.freqs, count.haps = count.haps, subst.model = subst.model, mu.rate = mu.rate, transi.rate = transi.rate, transv.rate = transv.rate)

HACSObj <- HACReal(p = 0.95, perms = 10000, conf.level = 0.95, num.iters = NULL, filename = "output")

### Run HAC simulations ###

set.seed(0673227) # set random seed for reproducibility, if desired
HAC.simrep(HACSObj) # simulate
