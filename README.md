# HACSim.R

A novel R simulation of haplotype accumulation curves. 

Species haplotypes are treated as distinct character labels ("1", "2", ...), where "1" denotes the most frequent haplotype, "2" denotes the second-most frequent haplotype, and so forth.

The present algorithm randomly samples species haplotype labels repeatedly in an iterative fashion, the idea being that levels of species haplotypic variation that are currently catalogued in the Barcode of Life Data Systems (BOLD; http://v4.boldsystems.org/) can serve as proxies for total haplotype diversity that may exist for a given species.

HACSim.R comprises two main functions (relevant code found in the files HAC.sim.R and HAC.simrep.R): 

> HAC.sim(K = 1, N, Hstar, probs, perms = 10000, p = 0.95, plot.out = FALSE)

> HAC.simrep().

Function arguments to HAC.sim() are as follows:

* **K** = Number of (equally-sized) (sub)populations/demes for a given species (**K** = 1 by default)

* **N** = Number of observed individuals (DNA sequences) of a given species 

* **Hstar** = Number of observed haplotypes (unique DNA sequences) for a given species

* **probs** = Haplotype frequency distribution for a given species

* **perms** = Number of permutations to generate species haplotype accumulation curve (**perms** = 10000 by default)

* **p** = Proportion of species haplotypes to recover (**p** = 0.95 by default)

* **plot.out** = An option to plot generated species haplotype accumulation curves (**plot.out** = FALSE by default)

**perms** controls the smoothness of generated haplotype accumulation curves. As **perms** &rarr; &infin;, haplotype accumulation curves "smooth out" and approach H* asymptotically.

HAC.sim() performs a single iteration of haplotype accumulation for a given species. If the desired level of haplotype recovery is not reached, then HAC.simrep() (which takes no arguments) is called in order to perform successive iterations until the desired fraction of haplotypes captured is at least **p**.

Both HAC.sim() and HAC.simrep() output simple "measures of closeness" for overall haplotype sampling completeness. Both absolute (counts) and relative (proportions) of species haplotypes sampled (observed) and missing (unobserved) are reported, along with estimates of the required sample size needed to uncover the specified level of species haplotypes and the number of additional specimens needed to be randomly sampled for a given species. In addition to haplotype accumulation curves, plots depicting empirical species haplotype frequency distributions are also displayed. 

Measures of closeness for overall sampling completeness are given by the following fomulae (Phillips *et al.*, 2015):

* Mean number of haplotypes sampled: H

* Number of haplotypes not sampled: H* - H

* Proportion of haplotypes sampled: H / H*

* Proportion of haplotypes not sampled: (H* - H) / H*

* Mean value of N*: NH* / H

* Mean number of individuals not sampled: N* - N

In addition to users specifying unique values for **N**, **Hstar** and **probs**, default parameters can also be altered in order to produce more interesting output (e.g., simulating multiple subpopulations). It may be necessary to increase **perms** in order to smooth out the curves, but this will increase algorithm runtime substantially. In addition, when there are multiple demes for a species (*i.e.*, **K** > 1), resulting statistical output will display sample size across *all* demes, while the graphical output will show haplotype accumulation curves and haplotype frequencies *per* deme. So, for example, if **N** = 50 and **K** = 2, then there would be **N**/**K** = 25 individuals in each deme.

To run the algorithm, do the following in a fresh R script:

1. Import required R scripts as follows:

> source("HAC.sim.R")

> source("HAC.simrep.R")

2. Set up all algorithm input parameters using standard R variable assignment. 

3. Run the following lines of code in succession:

> HAC.sim(K = K, N = N, Hstar = Hstar, probs = probs, r = r, perms = perms, p = p)

> HAC.simrep() **(needed only if HAC.sim() does not converge to at least p)**


### Examples (Step 2 above) ###

1. #### Equal haplotype frequency - Hypothetical species ####

* **N** = 100

* **Hstar** = 10

* **probs** = rep(1/Hstar, Hstar)

* **plot.out** = TRUE

2. #### Unequal haplotype frequency - Lake whitefish (*Coregonus clupeaformis*) #### see HAC.simrun.R

* **N** = 240
 
* **Hstar** = 15 
 
* **probs** = c(220/N, rep(3/N, 2), rep(2/N, 2), rep(1/N, 10))

* **plot.out** = TRUE

#### Custom user data ####

Users can implement their own custom species datasets mined from BOLD (not necessarily 5'-COI), but will first need to collapse DNA sequences into haplotypes and then extract the haplotype frequency distribution in order to determine values for **Hstar** and **probs**. This can be accomplished in many ways including, but not limited to the R package 'spider' (Brown *et al.*, 2012), or online interfaces such as FaBox (http://users-birc.au.dk/biopv/php/fabox/). In general, alignments should be of sufficiently high quality (*i.e.*, non-GenBank, and free of ambiguous/missing nucleotide bases, since these can lead to overestimation of overall haplotype diversity for a given species).  

Aligned and trimmed 652 bp COI barcode sequences for Lake whitefish (*Coregonus clupeaformis*) (Example 2 above) are included for download from this repository. 
