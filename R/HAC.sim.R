### Haplotype Accumulation Curve Simulation ###

##########

# Author: Jarrett D. Phillips
# Last modified: May 15, 2018

##########

## Best run in RStudio ##
## DO NOT change order of code (can throw errors)! ##

#####

## Input parameters ###

# Required #

# N = Number of specimens (DNA sequences)
# Hstar = Number of observed unique haplotypes
# probs = Probability frequency distribution of haplotypes
# perms = Number of permutations
# p = Proportion of unique haplotypes to recover

# Optional #

# input.seqs = Analyze inputted aligned/trimmed FASTA DNA sequence file (TRUE / FALSE)?
# prop.seqs = Proportion of DNA sequences to sample 



#####

HAC.sim <- function(N, 
                    Hstar, 
                    probs,
                    perms = 10000,
                    K = 1, # DO NOT CHANGE
                    p = 0.95,
                    input.seqs = FALSE,
                    subset.seqs = FALSE,
                    prop.seqs = NULL,
                    progress = TRUE) {
	
	cat("\n")
	
    if (progress == TRUE) {
      pb <- utils::txtProgressBar(min = 0, max = 1, style = 3)
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
		  
	  if (subset.seqs == TRUE) { # take random subset of sequences (e.g., prop.seqs = 0.10 (10%))
	                             # can be used to simulate migration/gene flow
	    seqs <- as.list(seqs)
		  seqs <- sample(seqs, size = ceiling(prop.seqs * length(seqs)), replace = FALSE)
		  seqs <- as.matrix(seqs)
	  }
		 
		assign("N", dim(seqs)[[1]], envir = .GlobalEnv)
		h <- sort(haplotype(seqs), decreasing = TRUE, what = "frequencies")
		rownames(h) <- 1:nrow(h)
		assign("Hstar", dim(h)[[1]], envir = .GlobalEnv)
		assign("probs", lengths(attr(h, "index")) / N, envir = .GlobalEnv)

	}
  
  ## Error messages ##
  
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
	
	num.specs <- N
		
	## Create an ID for each haplotype ##
	
	haps <- 1:Hstar
	
	## Assign individuals (N) ##
	
	specs <- 1:num.specs
	
	## Generate permutations. Assume each permutation has N individuals, and sample those 
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

  ## Update progress bar ##
	
    if (progress == TRUE) {
      utils::setTxtProgressBar(pb, i)
    }
	
	## Calculate the mean and CI for number of haplotypes recovered over all permutations

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
	  assign("X", ((N * Hstar) / P) - N, envir = .GlobalEnv)
	
	## Calculate slope of curve using last 10 points on curve
	
	# perms must be large enough to ensure monotonicity and a non-negative slope
    
    lin.reg <- lm(means ~ specs, data = tail(d, n = 10))
    b1 <- coef(lin.reg)[[2]]
	
    ## Output results to R console and text file ##
		    
	  cat("\n \n --- Measures of Sampling Closeness --- \n \n", 
	  "Mean number of haplotypes sampled: " , P, 
	  "\n Mean number of haplotypes not sampled: " , Q, 
	  "\n Proportion of haplotypes (specimens) sampled: " , R, 
	  "\n Proportion of haplotypes (specimens) not sampled: " , S,
	  "\n \n Mean value of N*: ", Nstar,
	  "\n Mean number of specimens not sampled: ", X, 
	  "\n \n Curve slope (last 10 points): ", b1,
	  "\n \n Mean number of specimens required to observe one new haplotype: ", 1 / b1, "\n \n")
	
    # name <- c("Mean number of haplotypes sampled", 
	         # "Mean number of haplotypes not sampled", 
	         # "Proportion of haplotypes (specimens) sampled", 
	         # "Proportion of haplotypes (specimens) not sampled",
	         # "Curve slope (last 10 points)", 
	         # "Mean number of specimens required to observe one new haplotype", 
	         # "Mean value of N*", 
	         # "Mean number of specimens not sampled")

  # write(P, file = "data.txt", append = TRUE)
  # write(Q, file = "Q.txt", append = TRUE)
  # write(R, file = "R.txt", append = TRUE)
  # write(S, file = "S.txt", append = TRUE)
  # write(b1, file = "b1.txt", append = TRUE)
  # write(1/b1, file = "invb1.txt", append = TRUE)
  # write(Nstar, file = "Nstar.txt", append = TRUE)
  # write(X, file = "X.txt", append = TRUE)
	# write.table(c(P, Q, R, S, b1, 1 / b1, Nstar, X), file = "data.txt", append = TRUE)
	
	## Plot the mean haplotype accumulation curve (averaged over perms number of curves) and haplotype frequency barplot ##

	par(mfrow = c(1, 2))
	plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes",  ylim = c(1, Hstar))
	polygon(x = c(specs, rev(specs)), y = c(lower, rev(upper)), col = "gray")
	lines(specs, means, lwd = 2)
	abline(v = N, lty = 2)
	abline(h = c(p * Hstar, R * Hstar), lty = 2)
	HAC.bar <- barplot(num.specs * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = haps)

}