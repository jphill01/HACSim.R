# HACSim #

[![Build Status](https://travis-ci.com/jphill01/HACSim.R.svg?token=R26mQLGR48w1Rk2BsVzG&branch=master)](https://travis-ci.com/jphill01/HACSim.R) [![Licence](https://img.shields.io/cran/l/HACSim.svg) [![](https://img.shields.io/cran/v/HACSim.svg) [![](https://cranlogs.r-pkg.org/badges/HACSim)](https://cran.r-project.org/package=HACSim)

HACSim (**H**aplotype **A**ccumulation **Sim**ulator) is a novel nonparametric stochastic (Monte Carlo) algorithm written in R for the simulation of haplotype accumulation curves. It can be employed to determine required sample sizes for DNA barcoding, specifically pertaining to recovery of total haplotype variation that may exist for a given species. 

Most DNA barcoding studies conducted to date suggest sampling between 5-10 individuals per species due to reseach costs. However, it has been shown that low sample sizes can greatly underestimate haplotype diversity for geograpically-widespread taxa. The present algorithm is in place to more accurately determine sample sizes that are needed to uncover all putative haplotypes that may exist for a given species. Implications of such an approach include accelerating the construction and growth of DNA barcode reference libraries for species of interest within the Barcode of Life Data Systems (BOLD; http://v4.boldsystems.org/).

Within the simulation algorithm, species haplotypes are treated as distinct character labels (1, 2, ...), where 1 denotes the most frequent haplotype, 2 denotes the second-most frequent haplotype, and so forth. The algorithm then randomly samples species haplotype labels in an iterative fashion, until all unique haplotypes have been observed. The idea is that levels of species haplotypic variation that are currently catalogued in BOLD can serve as proxies for total haplotype diversity that may exist for a given species.

Molecular loci besides DNA barcode genes (5'-COI, rbcL/matK, ITS regions) can be used with HACSim (*e.g.*, cyt*b*). As well, sequence databases other than BOLD can be targeted. 

HACSim comprises three main functions:

> HACHypothetical()

> HACReal()

> HAC.simrep().

Function arguments to HAC.sim() are passed from either HACHypothetical() or HACReal() as follows:

* **N** = Number of observed individuals (DNA sequences) of a given species 

* **Hstar** = Number of observed haplotypes (unique DNA sequences) for a given species

* **probs** = Haplotype frequency distribution vector for a given species (must have a length of **Hstar** and sum to 1)

* **perms** = Number of permutations to generate species haplotype accumulation curve (**perms** = 10000 by default)

* **p** = Proportion of species haplotypes to recover (**p** = 0.95 by default)

* **conf.level** = Desired confidence leval for graphical output and confidence interval calculation \\ (**conf.level** = 0.95 by default)

* **subsample** = Should a subsample of haplotype labels or DNA sequences be taken \\ (**subsample** = FALSE by default)

* **prop** = Proportion of haplotypw labels or DNA sequences to sample when **subsample** = TRUE \\ (**prop** = NULL by default)

* **progress** = Should iteration output be printed to the R console? (**progress** = TRUE by default)

* **filename** = Name of file where simulation results are to be saved (**filename** = NULL by default)

**perms** controls the smoothness of generated haplotype accumulation curves. As **perms** &rarr; &infin;, haplotype accumulation curves "smooth out" and approach H* asymptotically.

Resulting output for the first iteration of HAC.simrep() reflects current levels of sampling effort found within BOLD for a given species. If the desired level of haplotype recovery is not reached, then perform successive iterations are performed until the desired fraction of haplotypes captured is at least **p**.

Setting **p** = 0.95 corresponds to uncovering 95% of all haplotypes that may exist for a given species. At this level, the generated haplotype accumulation curve reaches a slope close to zero and further sampling effort is unlikely to uncover any new haplotypes. 

HAC.simrep() outputs simple "Measures of Sampling Closeness" for overall haplotype sampling completeness. Both absolute (counts) and relative (proportions) of species haplotypes sampled (observed) and missing (unobserved) are reported, along with estimates of the required sample size needed to uncover the specified number of species haplotypes and the number of additional specimens needed to be randomly sampled for a given species. In addition to haplotype accumulation curves, plots depicting species haplotype frequency distributions are also displayed. 

Measures of Sampling Closeness for overall haplotype/specimen sampling completeness are given by the following fomulae (Phillips *et al.*, 2015):

* Mean number of haplotypes sampled: H 

* Number of haplotypes not sampled: H* - H

* Proportion of haplotypes (specimens) sampled: H / H*

* Proportion of haplotypes (specimens) not sampled: (H* - H) / H* 

* Mean value of N*: NH* / H

* Mean number of specimens not sampled: N* - N

In addition to users specifying unique values for **N**, **Hstar** and **probs**, default parameters can also be altered in order to produce more interesting output (e.g., subsampling of haplotype labels or DNA sequences). It may be necessary to increase **perms** in order to smooth out the curves, but this will increase algorithm runtime substantially. 

### Running the Simulation ###

To run the algorithm, do the following in a fresh R script:

1. Set up all desired algorithm input parameters using standard R variable assignment. 

2. Run either of the following lines of code for simulation of either hypothetical or real species:
> HACHypothetical(...) 

> HACReal(...)

3. Run the simulator
> HAC.simrep(...)

**NOTE**: Users must have a system compiler (e.g. Xcode Command Line Tools on MacOS) installed in order to run the algorithm successfully.  

Depending on the size of input parameters to the simulation, the algorithm in its current form, can be quite slow to reach full convergence. For a species with equal haplotype frequency, saturation of the haplotype accumulation curve is (usually) reached very rapidly (one exception occurs when **N** = **Hstar**). On the other hand, for a species with many rare haplotypes, the generated haplotype accumulation curve will take signficantly longer to reach an asymptote, since rare haplotypes will not be sampled as frequently as dominant ones.

To sucessfully run the simulation algorithm, the following conditions must hold:

* **N** must be greater than 1

* **N** must be greater than or equal to **Hstar**

* **probs** must sum to 1 and must have a length of **Hstar**

### Examples (Step 2 above) ###

1. #### Equal haplotype frequency - Hypothetical species ####

* **N** = 100

* **Hstar** = 10

* **probs** = rep(1/Hstar, Hstar)

> HACHypothetical(N, Hstar, probs, ...)

2. #### Unequal haplotype frequency - Lake whitefish (*Coregonus clupeaformis*) ####

* **N** = 235
 
* **Hstar** = 15 
 
* **probs** = c(215/**N**, rep(3/**N**, 2), rep(2/**N**, 2), rep(1/**N**, 10)) (or see **Custom User Data** below)

> HACReal()

### Custom User Data ###

Users can implement their own custom species barcode datasets mined from BOLD (not necessarily 5'-COI) through setting the argument **input.seqs** to TRUE. Values for **N**, **Hstar** and **probs** are calculated automatically via the R packages 'ape' (Analyses of Phylogenetics and Evolution; Paradis, 2004) and 'pegas' (Phylogenetics and Evolution in R; Paradis *et al.* (2010)), so there is no need to specify these parameters before running HAC.simrep(). When HAC.simrep() is run, a pop-up window will appear prompting the user to select a previously aligned/trimmed FASTA DNA sequence file.

**NOTE**: Inputted sequence alignment files should be free of missing/ambiguous nucleotide data (Ns, gaps (--) and IUPAC ambiguity codes) in order to avoid overestimation of intraspecific haplotype diversity. 

Users can also subsample a random subset of DNA sequences or haplotypes with the arguments **subsample** and **prop**. This can be employed to simulate migration/gene flow of individuals/haplotypes, or to reduce computation overhead. 

Aligned and trimmed 652 bp 5'-COI barcode sequences for Lake whitefish (*Coregonus clupeaformis*) (**Example 2** above) are included for download from this repository.

### Algorithm Availability ###

HACSim is available for download as a package from the Comprehensive R Archive Network (CRAN).

> install.packages(HACSim)

> library(HACSim).

The reference manual for HACSim, which includes built-in functions with explanations for their proper use, can be accessed by typing
> ?HACSim.

### Citing This Work ###

In utilizing HACSim, the following publications may be of interest and worth citing:

**Phillips, J.D.**, Gwiazdowski, R.A., Ashlock, D. and Hanner, R. (2015). An exploration of sufficient sampling effort to describe intraspecific DNA barcode haplotype diversity: examples from the ray-finned fishes (Chordata: Actinopterygii). *DNA Barcodes*, **3**: 66-73.

**Phillips, J.D.**, Gillis, D.J. and Hanner, R.H. (2019). Incomplete estimates of genetic diversity within species: Implications for DNA barcoding. *Ecology and Evolution*, 1-15.

