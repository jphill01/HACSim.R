# HACSim.R

R simulation of haplotype accumulation curves.

HACSim.R comprises two main functions: HAC.sim(K = 1, N, Hstar, probs, m = 0, perms = 10000, p = 0.95) and HAC.simrep().

Function arguments to HAC.sim() are as follows:

K = number of (sub)populations/demes (K = 1 by default)
N = Number of obsered individuals (DNA barcode sequences) of a given species 
Hstar = Number of observed haplotypes (unique DNA barcode sequences) for a given species
probs = Haplotype probability distribution for a given species
m = Overall migration rate of individuals/haplotypes between demes (m = 0 by default)
perms = Number of permutations to generate species haplotype accumulation curve (perms = 10000 by default)
p = Proportion of species haplotypes to recover (p = 0.95 by default)
