# HACSim.R

A novel R simulation of haplotype accumulation curves.

HACSim.R comprises two main functions: HAC.sim(K = 1, N, Hstar, probs, m = 0, perms = 10000, p = 0.95) and HAC.simrep().

HAC.sim() performs a single iteration of haplotype accumulation for a given species. If the desired level of haplotype recovery is not reached, then HAC.simrep() (which takes no arguments) is called in order to perforn successive iterations until the desired fraction of haplotypes captured is at least p.


Function arguments to HAC.sim() are as follows:

K = Number of (sub)populations/demes for a given species (K = 1 by default)

N = Number of observed individuals (DNA sequences) of a given species 

Hstar = Number of observed haplotypes (unique DNA sequences) for a given species

probs = Haplotype frequency distribution for a given species

m = Overall migration rate of individuals/haplotypes between demes (m = 0 by default)

perms = Number of permutations to generate species haplotype accumulation curve (perms = 10000 by default)

p = Proportion of species haplotypes to recover (p = 0.95 by default)

Both HAC.sim() and HAC.simrep() output simple "measaures of closeness" for overall haplotype sampling completeness. Both absolute (counts) and relative (proprtions) of species haplotypes sampled (observed) and missing (unobserved) are reported, along with estimates of the required sample size needed to uncover the specified level of species haplotypes and the number of additional specimens needed to be randomly sampled for a given species 

In addition to users specifying unique values for N, Hstar and probs, default parameters can also be altered in order to produce more interesting output. 

# Example: Lake Whitefish ($Coregonus clupeaformis$) #

