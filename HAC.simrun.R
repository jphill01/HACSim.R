##### Unequal haplotype frequency - Lake whitefish (Coregonus clupeaformis) #####

### Set parameters ###

K <- 1  # number of equally-sized (sub)populations
N <- 240 # total number of sampled individuals
Hstar <- 15 # total number of haplotypes
probs <- c(220/N, rep(3/N, 2), rep(2/N, 2), rep(1/N, 10)) # haplotype frequency distribution
perms <- 10000 # number of permutations
p <- 0.95 # proportion of haplotypes to recover 
plot.out <- TRUE # plot curve and barplot


### Run simulations ###

ptm <- proc.time() # set timer

HAC.sim(K = K,  N = N, Hstar = Hstar, probs = probs, perms = perms, p = p, plot.out = plot.out)
HAC.simrep()

proc.time() - ptm
