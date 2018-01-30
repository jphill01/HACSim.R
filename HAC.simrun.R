##### Unequal haplotype frequency - Lake whitefish (Coregonus clupeaformis) #####

### Set working directory ###

setwd("/Users/jarrettphillips/Desktop/HAC simulation")

### Clear memory ###

remove(list = ls())

### Load libraries ###

library(Rcpp)
library(RcppArmadillo)

### Load scripts ###

source("HACaccum.cpp")
source("migrate.cpp")
source("HAC.sim.R")
source("HAC.simrep.R")

### Set parameters ###

N <- 240 # total number of sampled individuals
Hstar <- 15 # total number of haplotypes
probs <- c(220/N, rep(3/N, 2), rep(2/N, 2), rep(1/N, 10)) # haplotype frequency distribution
K <- 1  # number of equally-sized (sub)populations
m <- 0 # migration rate between subpopulations
model <- NULL
perms <- 10000 # number of permutations
p <- 1 # proportion of haplotypes to recover 
plot.out <- TRUE # haplotype accumulation curve and haplotype frequency barplot


### Run simulations ###

ptm <- proc.time() # set timer

HAC.sim(N = N, Hstar = Hstar, probs = probs, K = K, m = m, model = model, perms = perms, p = p, plot.out = plot.out)
HAC.simrep()

proc.time() - ptm
