# HACSim.R #

HACSim.R is a novel R simulation of haplotype accumulation curves. The present simulation algorithm can be employed in order to determine required sample sizes for DNA barcoding, specifically as pertaining to recovery of total haplotype variation for a given species.

Most DNA barcoding studies conducted to date suggest sampling between 5-10 individuals per species due to reseach costs. However, it has been shown that low sample sizes can greatly underestimate haplotype diversity for geograpically-widespread taxa. The present algorithm is in place to more accurately determine sample sizes that are needed to uncover all putative haplotypes that may exist for a given species. Implications of such an approach include accelerating the construction and growth of DNA barcode reference libraries for species of interest within the Barcode of Life Data Systems (BOLD; http://v4.boldsystems.org/).

Within the simulation algorithm, species haplotypes are treated as distinct character labels (1, 2, ...), where 1 denotes the most frequent haplotype, 2 denotes the second-most frequent haplotype, and so forth. The algorithm then randomly samples species haplotype labels in an iterative fashion, until all unique haplotypes have been observed at least once. The idea is that levels of species haplotypic variation that are currently catalogued in BOLD can serve as proxies for total haplotype diversity that may exist for a given species.

HACSim.R comprises two main functions (relevant code found in the files HAC.sim.R and HAC.simrep.R):

> HAC.sim()

> HAC.simrep().

Function arguments to HAC.sim() are as follows:

* **N** = Number of observed individuals (DNA sequences) of a given species 

* **Hstar** = Number of observed haplotypes (unique DNA sequences) for a given species

* **probs** = Haplotype frequency distribution for a given species (must sum to 1)

* **perms** = Number of permutations to generate species haplotype accumulation curve (**perms** = 10000 by default)

* **p** = Proportion of species haplotypes to recover (**p** = 0.95 by default)

* **input.seqs** = Analyze a user-specified DNA sequence FASTA file?
(**input.seqs** = FALSE by default)

* **subset.seqs** = Subsample DNA sequences? (**subset.seqs** = FALSE by default)

* **subset.haps** = Subsample DNA sequences? (**subset.haps** =  NULL by default)

* **prop.seqs** = Proportion of DNA sequences to subsample (**prop.seqs** = NULL by default)

* **prop.haps** = Proportion of haplotypes to subsample (**prop.haps** = NULL by default)

* **num.pts** = Number of terminal data points used to calculate curve slope (**num.pts = 10** by default)

* **prop.pts** = Proportion of terminal data points used to compute curve slope (**prop.pts** = NULL by default)

**perms** controls the smoothness of generated haplotype accumulation curves. As **perms** &rarr; &infin;, haplotype accumulation curves "smooth out" and approach H* asymptotically.

HAC.sim() performs a single iteration of haplotype accumulation for a given species. Resulting output reflects current levels of sampling effort found within BOLD for a given species. If the desired level of haplotype recovery is not reached, then HAC.simrep() (which takes no arguments) is called in order to perform successive iterations until the desired fraction of haplotypes captured is at least **p**.

Setting **p** = 0.95 corresponds to uncovering 95% of all haplotypes that may exist for a given species. At this level, the generated haplotype accumulation curve reaches a slope close to zero and further sampling effort is unlikely to uncover any new haplotypes. 

Users can specify either **num.pts** or **prop.pts** in calculating the terminal slope of the generated haplotype accumulation curve.

Both HAC.sim() and HAC.simrep() output simple "Measures of Sampling Closeness" for overall haplotype sampling completeness. Both absolute (counts) and relative (proportions) of species haplotypes sampled (observed) and missing (unobserved) are reported, along with estimates of the required sample size needed to uncover the specified level of species haplotypes and the number of additional specimens needed to be randomly sampled for a given species. In addition to haplotype accumulation curves, plots depicting species haplotype frequency distributions are also displayed. 

Measures of Sampling Closeness for overall haplotype/specimen sampling completeness are given by the following fomulae (Phillips *et al.*, 2015):

* Mean number of haplotypes sampled: H

* Number of haplotypes not sampled: H* - H

* Proportion of haplotypes (specimens) sampled: H / H* (alternatively, N / N*)

* Proportion of haplotypes (specimens) not sampled: (H* - H) / H* (alternatively, (N* - N) / N*)

* Mean value of N*: NH* / H

* Mean number of specimens not sampled: N* - N

Generated output also returns the slope estimate (*&beta;*<sub>1</sub>) for the haplotype accumulation curve, calculated using the last 10 points occurring on the curve. In addition, the quantity 1 / *&beta;*<sub>1</sub> is also returned. This is a measure of the number of additional specimens needed to be sampled to observe one new haplotype.

In addition to users specifying unique values for **N**, **Hstar** and **probs**, default parameters can also be altered in order to produce more interesting output (e.g., simulating multiple subpopulations with or without migration/gene flow). It may be necessary to increase **perms** in order to smooth out the curves, but this will increase algorithm runtime substantially. 

### Running the Simulation ###

To run the algorithm, do the following in a fresh R script:

1. Import required C++/R scripts as follows:

> source("accumulate.cpp")

> source("HAC.sim.R")

> source("HAC.simrep.R")

2. Set up all algorithm input parameters using standard R variable assignment. 

3. Run the following line of code:

> HAC.simrep()

**NOTE**: Users must have a system compiler (e.g. Xcode Command Line Tools on MacOS) installed in order to run the algorithm sucessfully.  

Depending on the size of input parameters to the simulation, the algorithm in its current form, can be quite slow to reach full convergence. For a species with equal haplotype frequency, saturation of the haplotype accumulation curve is (usually) reached very rapidly (one exception occurs when **N** = **Hstar**). On the other hand, for a species with many rare haplotypes, the generated haplotype accumulation curve will take signficantly longer to reach an asymptote, since rare haplotypes will not be sampled as frequently as dominant ones.

In order to sucessfully run the simulation algorithm, the following conditions must hold:

* **N** must be greater than 1

* **N** must be greater than or equal to **Hstar**

* **probs** must sum to 1

### Examples (Step 2 above) ###

1. #### Equal haplotype frequency - Hypothetical species ####

* **N** = 100

* **Hstar** = 10

* **probs** = rep(1/Hstar, Hstar)

2. #### Unequal haplotype frequency - Lake whitefish (*Coregonus clupeaformis*) ####

* **N** = 240
 
* **Hstar** = 15 
 
* **probs** = c(214/**N**, rep(3/**N**, 2), rep(2/**N**, 2), rep(1/**N**, 10)) (or see **Custom User Data** below)


### Custom User Data ###

Users can implement their own custom species barcode datasets mined from BOLD (not necessarily 5'-COI) through setting the argument **input.seqs** to TRUE. Values for **N**, **Hstar** and **probs** are calculated automatically via the R packages 'ape' (Analyses of Phylogenetics and Evolution; Paradis, 2004) and 'pegas' (Phylogenetics and Evolution in R; Paradis *et al.* (2010)), so there is no need to specify these parameters before running HAC.simrep(). When HAC.simrep() is run, a pop-up window will appear prompting the user to select a previously aligned/trimmed FASTA DNA sequence file.

**NOTE**: Inputted sequence alignment files should be free of missing/ambiguous nucleotide data (Ns, gaps (--) and IUPAC ambiguity codes) in order to avoid overestimation of intraspecific haplotype diversity. 

Users can also subsample a random subset of DNA sequences or haplotypes with the arguments (**subset.seqs**, **prop.seqs** **subset.haps**, and **prop.haps**). This can be employed to simulate migration/gene flow of individuals/haplotypes. 

Aligned and trimmed 652 bp 5'-COI barcode sequences for Lake whitefish (*Coregonus clupeaformis*) (**Example 2** above) are included for download from this repository.

### Citing This Work ###

In utilizing HACSim, please cite the following publication:

**Phillips, J.D.**, Gwiazdowski, R.A., Ashlock, D. and Hanner, R. (2015). An exploration of sufficient sampling effort to describe intraspecific DNA barcode haplotype diversity: examples from the ray-finned fishes (Chordata: Actinopterygii). *DNA Barcodes*, **3**: 66-73. 
