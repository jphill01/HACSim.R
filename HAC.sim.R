### Haplotype Accumulation Curve Simulation ###

##########

# Author: Jarrett D. Phillips
# Last modified: April 2, 2018

##########

## Best run in RStudio ##
## DO NOT change order of code (can throw errors)! ##

#####

## Input parameters ###

# N = Number of specimens (DNA sequences)
# Hstar = Number of observed unique haplotypes
# probs = Probability frequency distribution of haplotypes
# K = Number of (sub)populations (demes, sampling sites)
# m = Overall migration rate between demes
# perms = Number of permutations
# p = Proportion of unique haplotypes to recover
# input.seqs = Analyze inputted aligned/trimmed FASTA DNA sequence file (TRUE / FALSE)?


#####

HAC.sim <- function(N, 
                    Hstar, 
                    probs, 
                    K = 1, 
                    m = 0,
                    perms = 10000, 
                    p = 0.95, 
                    input.seqs = FALSE, 
                    progress = TRUE) {
	
	cat("\n")
	
	if (progress == TRUE) {
		pb <- utils::txtProgressBar(min = 0, max = K, style = 3)
	}

	## Load DNA sequence data and set N, Hstar and probs ##
	
	if (input.seqs == TRUE) {
		seqs <- read.dna(file = file.choose(), format = "fasta")
		if (all(base.freq(seqs, all = TRUE)[5:17] != 0)) {
			warning("Inputted DNA sequences contain missing and/or ambiguous nucleotides, which may lead 
			        to overestimation of the number of observed unique haplotypes.  Consider excluding 
			        sequences or alignment sites containing these data. If missing and/or ambiguous bases 
			        occur at the ends of sequences, further alignment trimming is an option.")
		}
		assign("N", dim(seqs)[[1]], envir = .GlobalEnv)
		h <- sort(haplotype(seqs), decreasing = TRUE, what = "frequencies")
		rownames(h) <- 1:nrow(h)
		assign("Hstar", dim(h)[[1]], envir = .GlobalEnv)
		assign("probs", lengths(attr(h, "index")) / N, envir = .GlobalEnv)

	}
  
  ## Error messages ##
  
  if (N < K) {
    stop("N must be greater than or equal to K")
  }
  
  if (N < Hstar) {
    stop("N must be greater than or equal to Hstar")
  }
  
  if (N == 1) {
    stop("N must be greater than 1")
  }
  
  if (sum(probs) != 1) {
    stop("probs must sum to 1")
  }
	
	## Set up container(s) to hold the identity of each individual from each permutation ##
	
	num.specs <- ceiling(N / K)
		
	## Create an ID for each haplotype ##
	
	haps <- 1:Hstar
	
	## Assign individuals (N) to each subpopulation (K) ##
	
	specs <- 1:num.specs
	
	## Generate permutations, assume each permutation has N/K individuals, and sample those 
	# individuals' haplotypes from the probabilities ##
	
	gen.perms <- function() {
		sample(haps, size = num.specs, replace = TRUE, prob = probs)
	}
	
	pop <- array(dim = c(perms, num.specs, K))
	
	for (i in 1:K) {
		pop[,, i] <- replicate(perms, gen.perms())
	}

	## Perform haplotype accumulation ##
	
	HAC.mat <- accumulate(pop, specs, perms, K)
	
	## Allow individuals to migrate between subpopulations according to migration rate m ##
	
	if (K > 1 && m != 0) {
	  
	  inds1 <- sample(perms, size = ceiling(num.specs * m), replace = TRUE)
	  inds2 <- sample(perms, size = ceiling(num.specs * m), replace = TRUE)
	  
	  for (i in 1:K) {
	    for(j in 1:K) {
	      tmp <- pop[inds1,, i]
	      pop[inds1,, i] <- pop[inds2,, j]
	      pop[inds2,, j] <- tmp
	    }
	  }
	  
	}
	
	## Update progress bar ##
	
	if (progress == TRUE) {
	  utils::setTxtProgressBar(pb, i)
	}
	
	## Calculate the mean and CI for number of haplotypes recovered

	means <- apply(HAC.mat, MARGIN = 2, mean)
	lower <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.025))
	upper <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.975))
	
	## Make data accessible to user ##
	
	assign("d", data.frame(specs, means), envir = .GlobalEnv)
	
	## Compute simple summary statistics and display output ##
	
	# tail() is used here instead of max() because curves will not be monotonic if perms is not set high enough. When perms is large (say 10000), tail() is sufficiently close to max()  
	
	P <- tail(means, n = 1)
	Q <- Hstar - P
	assign("R", P / Hstar, envir = .GlobalEnv)
	S <- (Hstar - P) / Hstar
	assign("Nstar", (N * Hstar) / P, envir = .GlobalEnv)
	X <- ((N * Hstar) / P) - N
	
	## Calculate slope of curve using last 10 points on the curve ##
	
	# perms must be large enough to ensure monotonicity and a non-negative slope
	
	lin.reg <- lm(means ~ specs, data = tail(d, n = 10))
	b1 <- coef(lin.reg)[[2]]
	
	# Compute Nei's Haplotype diversity
	
	hd <- (Nstar / (Nstar - 1)) * (1 - sum(probs^2))
	
	## Output results to R console and text file ##
		    
	cat("\n \n --- Measures of Sampling Closeness --- \n \n", 
	"Mean number of haplotypes sampled: " , P, 
	"\n Mean number of haplotypes not sampled: " , Q, 
	"\n Proportion of haplotypes (specimens) sampled: " , R, 
	"\n Proportion of haplotypes (specimens) not sampled:  " , S, 
	"\n \n Haplotype diversity: ", hd, 
	"\n \n Curve slope (last 10 points): ", b1,
	"\n \n Mean number of specimens required to observe one new haplotype: ", 1 / b1,
	"\n \n Mean value of N*: ", Nstar / K, 
	"\n Mean number of specimens not sampled: ", X / K, "\n \n")
	
	name <-c("Mean number of haplotypes sampled", "Mean number of haplotypes not sampled", 
	         "Proportion of haplotypes (specimens) sampled", "Proportion of haplotypes (specimens) not sampled", 
	         "Haplotype diversity", "Curve slope (last 10 points)", "Mean number of specimens required to observe one new haplotype", 
	         "Mean value of N*", "Mean number of specimens not sampled")
	tbl <- c(P, Q, R, S, Nstar / K, X / K, b1, 1 / b1, hd)
	write(tbl, "data.txt", sep = "\n", append = TRUE)
	
	## Plot the haplotype accumulation curve and haplotype frequency barplot ##

			par(mfrow = c(1, 2))
			plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes",  ylim = c(1, Hstar))
			polygon(x = c(specs, rev(specs)), y = c(lower, rev(upper)), col = "gray")
			lines(specs, means, lwd = 2)
			HAC.bar <- barplot(num.specs * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = 1:Hstar)
			
}
