### Haplotype Accumulation Curve Simulation ###

##########

# Author: Jarrett D. Phillips
# Last modified: November 13, 2018

##########

## Best run in RStudio ##
## DO NOT change order of code (can throw errors)! ##

#####

## Input parameters ###

# Required #

# N = Number of specimens (DNA sequences)
# Hstar = Number of observed unique haplotypes
# probs = Probability frequency distribution of haplotypes

# Optional #

# p = Proportion of unique haplotypes to recover
# perms = Number of permutations (replications)
# input.seqs = Analyze inputted aligned/trimmed FASTA DNA sequence file (TRUE / FALSE)?
# sim.seqs = Simulate DNA sequences (TRUE / FALSE)?
# subst.model = Nucleotide substition model (JC69 / K80 / F81 / HKY85)
# subset.seqs = Subset of DNA sequences to sample
# prop.seqs = Proportion of DNA sequences to sample 
# prop.haps = Proportion of haplotypes to sample 
# subset.haps = Subset of haplotypes to sample
# num.pts = Number of points used to calculate curve slope 
# prop.pts = Proportion of points used to calculate curve slope

#####

HAC.sim <- function(N, 
                    Hstar, 
                    probs,
                    perms = 10000,
                    K = 1, # DO NOT CHANGE,
                    p = 0.95,
                    subset.haps = NULL,
                    prop.haps = NULL,
                    input.seqs = FALSE,
                    sim.seqs = FALSE,
                    subst.model = NULL,
                    subset.seqs = FALSE,
                    prop.seqs = NULL,
                    df = NULL, # dataframe
                    progress = TRUE) {
  
    cat("\n \n")

  ## Display progress bar ##
    
    if (progress == TRUE) {
      pb <- utils::txtProgressBar(min = 0, max = 1, style = 3)
    }

	## Load DNA sequence data and set N, Hstar and probs ##
	  
    if ((input.seqs == TRUE) || (sim.seqs == TRUE)) {
		  seqs <- read.dna(file = file.choose(), format = "fasta")
		  
		  bf <- base.freq(seqs, all = TRUE)[5:17]
	  
		if (any(bf > 0)) {
      stop("Inputted DNA sequences contain missing and/or ambiguous 
	    nucleotides, which may lead to overestimation of the number of 
	    observed unique haplotypes. Consider excluding sequences or alignment 
	    sites containing these data. If missing and/or ambiguous bases occur 
	    at the ends of sequences, further alignment trimming is an option.")
		}
		
		ptm <<- proc.time()
		  
	  if (subset.seqs == TRUE) { # take random subset of sequences (e.g., prop.seqs = 0.10 (10%))
	                             # can be used to simulate migration/gene flow
		  seqs <- seqs[sample(nrow(seqs), size = ceiling(prop.seqs * nrow(seqs)), replace = FALSE), ]
		  write.dna(seqs, file = "seqs.fas", format = "fasta")
	  }
		 
		assign("N", dim(seqs)[[1]], envir = .GlobalEnv)
		h <- sort(haplotype(seqs), decreasing = TRUE, what = "frequencies")
		rownames(h) <- 1:nrow(h)
		assign("Hstar", dim(h)[[1]], envir = .GlobalEnv)
		assign("probs", lengths(attr(h, "index")) / N, envir = .GlobalEnv)

	  }

  ## Simulate DNA sequences ##
  
    if (sim.seqs == TRUE) {
      
      dis <- dist.dna(seqs, model = subst.model) # K2P divergences
      
      tr <- nj(dis) # unrooted neighbour-joining tree
      
      bf <- base.freq(seqs)
      
      res <- as.DNAbin(simSeq(tr, l = ncol(seqs), bf = bf)) # convert to DNAbin format
    
      write.dna(res, file = "res.fas", format = "fasta") # output sequences to file
    
      assign("N", dim(res)[[1]], envir = .GlobalEnv)
      h <- sort(haplotype(res), decreasing = TRUE, what = "frequencies")
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
  
  ## Set up container to hold the identity of each individual from each permutation ##
    
    num.specs <- N
		
	## Create an ID for each haplotype ##
	  
    if (is.null(subset.haps)) {
	    haps <- 1:Hstar
	  } else {
	    subset.haps <- subset.haps
	  }
	  
	## Assign individuals (N) ##
	  
    specs <- 1:num.specs
	
	## Generate permutations. Assume each permutation has N individuals, and sample those 
	## individuals' haplotypes from the probabilities ##
	  
    gen.perms <- function() {
	    if (is.null(subset.haps)) {
	      sample(haps, size = num.specs, replace = TRUE, prob = probs)
	    } else {
	      resample <- function(x, ...) x[sample.int(length(x), ...)]
	      resample(subset.haps, size = num.specs, replace = TRUE, prob = probs[subset.haps])
	    }
	  }
	  
	  pop <- array(dim = c(perms, num.specs, K))
	  
	  for (i in 1:K) {
	    pop[,, i] <- replicate(perms, gen.perms())
	  }
	  
	## Perform haplotype accumulation ##
    
	  HAC.mat <- accumulate(pop, specs, perms, K) # one row selected at random 

  ## Update progress bar ##
    
	  if (progress == TRUE) {
      utils::setTxtProgressBar(pb, i)
    }
	
	## Calculate the mean and CI for number of haplotypes recovered over all permutations
	  
	  means <- apply(HAC.mat, MARGIN = 2, mean)
	  lower <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.025)) 
	  upper <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.975)) 
	
	## Make data accessible to user ##
	 
	  assign("d", data.frame(specs, means, lower, upper), envir = .GlobalEnv)

	## Compute simple summary statistics and display output ##
	## tail() is used here instead of max() because curves will not be monotonic if perms is not set high enough. When perms is large (say 10000), tail() is sufficiently close to max()
	 
	  P <- tail(means, n = 1)
	
	 if (is.null(subset.haps)) {
	   Q <- Hstar - P
	   assign("R", P / Hstar, envir = .GlobalEnv)
	   S <- (Hstar - P) / Hstar
	   assign("Nstar", (N * Hstar) / P, envir = .GlobalEnv)
	   X <- ((N * Hstar) / P) - N
	 } else {
	   Q <- length(subset.haps) - P
	   assign("R", P / length(subset.haps), envir = .GlobalEnv)
	   S <- (length(subset.haps) - P) / length(subset.haps)
	   assign("Nstar", (N * length(subset.haps)) / P, envir = .GlobalEnv)
	   X <- ((N * length(subset.haps)) / P) - N
	 }
	
  ## Output results to R console and CSV file ##
	   
	   cat("\n \n --- Measures of Sampling Closeness --- \n \n", 
	       "Mean number of haplotypes sampled: " , P, "( 95% CI:", paste(ceiling(max(lower)), ceiling(max(upper)), sep = "-"), ")",
	       "\n Mean number of haplotypes not sampled: " , Q, 
	       "\n Proportion of haplotypes sampled: " , R, 
	       "\n Proportion of haplotypes not sampled: " , S,
	       "\n \n Mean value of N*: ", Nstar,
	       "\n Mean number of specimens not sampled: ", X)

    df[nrow(df) + 1, ] <- c(P, ceiling(max(lower)), ceiling(max(upper)), Q, R, S, Nstar, X)
    
  ## Plot the mean haplotype accumulation curve (averaged over perms number of curves) and haplotype frequency barplot ##
      
    par(mfrow = c(1, 2))
    
      if (is.null(subset.haps)) {
        plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes",  ylim = c(1, Hstar), main = "Haplotype accumulation curve")
      } else {
        plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes",  ylim = c(1, length(subset.haps)), main = "Haplotype accumulation curve")
      }
      
      polygon(x = c(specs, rev(specs)), y = c(lower, rev(upper)), col = "gray")
      lines(specs, means, lwd = 2)
      
      if (is.null(subset.haps)) {
        abline(h = R * Hstar, v = N, lty = 2) # dashed line
        abline(h = p * Hstar, lty = 3) # dotted line
        HAC.bar <- barplot(num.specs * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = haps, main = "Haplotype frequency distribution")
      } else {
        abline(h = R * length(subset.haps), v = N, lty = 2) # dashed line
        abline(h = p * length(subset.haps), lty = 3) # dotted line
        HAC.bar <- barplot(num.specs * (probs[subset.haps] / sum(probs[subset.haps])), xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = subset.haps, main = "Haplotype frequency distribution")
      }
      
	  df
} # end HAC.sim
