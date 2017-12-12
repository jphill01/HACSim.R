# HACSim.R

A novel R simulation of haplotype accumulation curves.

Species haplotypes are treated as distinct character labels ("1", "2", ...), where "1" denotes the most frequent haplotype, "2" denotes the second-most frequent haplotype, and so forth.

HACSim.R comprises two main functions: HAC.sim(K = 1, N, Hstar, probs, m = 0, perms = 10000, p = 0.95) and HAC.simrep().

HAC.sim() performs a single iteration of haplotype accumulation for a given species. If the desired level of haplotype recovery is not reached, then HAC.simrep() (which takes no arguments) is called in order to perform successive iterations until the desired fraction of haplotypes captured is at least p.

Function arguments to HAC.sim() are as follows:

* K = Number of (sub)populations/demes for a given species (K = 1 by default)

* N = Number of observed individuals (DNA sequences) of a given species 

* Hstar = Number of observed haplotypes (unique DNA sequences) for a given species

* probs = Haplotype frequency distribution for a given species

* m = Overall migration rate of individuals/haplotypes between demes (m = 0 by default)

* perms = Number of permutations to generate species haplotype accumulation curve (perms = 10000 by default)

* p = Proportion of species haplotypes to recover (p = 0.95 by default)

perms controls the smoothness of generated haplotype accumulation curves. 

Both HAC.sim() and HAC.simrep() output simple "measaures of closeness" for overall haplotype sampling completeness. Both absolute (counts) and relative (proprtions) of species haplotypes sampled (observed) and missing (unobserved) are reported, along with estimates of the required sample size needed to uncover the specified level of species haplotypes and the number of additional specimens needed to be randomly sampled for a given species 

In addition to users specifying unique values for N, Hstar and probs, default parameters can also be altered in order to produce more interesting output. perms should not be altered 

### Example ###

Lake whitefish (*Coregonus clupeaformis*)

 N = 240
 
 Hstar = 15
 
 probs = c(220/N, rep(3/N, 2), rep(2/N, 2), rep(1/N, 10))
