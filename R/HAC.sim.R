### Haplotype Accumulation Curve Simulation ###

##########

# Author: Jarrett D. Phillips
# Last modified: August 2, 2018

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
                    num.seqs = NULL,
                    length.seqs = NULL,
                    transition.rate = NULL,
                    transversion.rate = NULL,
                    subset.seqs = FALSE,
                    prop.seqs = NULL,
                    num.pts = 10,
                    prop.pts = NULL,
                    df = NULL, # dataframe
                    progress = TRUE) {
	
	  cat("\n")
  
  ## Display progress bar ##
	
    if (progress == TRUE) {
      pb <- utils::txtProgressBar(min = 0, max = K, style = 3)
    }

	## Load DNA sequence data and set N, Hstar and probs ##
	
	  if (input.seqs == TRUE) {
		  seqs <- read.dna(file = file.choose(), format = "fasta")
	  if (any(base.freq(seqs, all = TRUE)[5:17] > 0)) {
		  stop("Inputted DNA sequences contain missing and/or ambiguous 
	nucleotides, which may lead to overestimation of the number of 
	observed unique haplotypes. Consider excluding sequences or alignment 
	sites containing these data. If missing and/or ambiguous bases occur 
	at the ends of sequences, further alignment trimming is an option.")
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
  
  if (sim.seqs == TRUE) {
    
    NUCL <- as.DNAbin(c("a","t","c","g"))
    TRANSISET <- list('a'=as.DNAbin('g'), 'g'=as.DNAbin('a'), 'c'=as.DNAbin('t'), 't'=as.DNAbin('c'))
    TRANSVSET <- list('a'=as.DNAbin(c('c','t')), 'g'=as.DNAbin(c('c','t')), 'c'=as.DNAbin(c('a','g')), 't'=as.DNAbin(c('a','g')))
    
    ## AUXILIARY FUNCTIONS ##
    ## generate sequence from scratch
    
    res <- sample(NUCL, size=length.seqs, replace=TRUE)
    
    ## create transitions for defined SNPs
    transi <- function(snp){
      res <- unlist(TRANSISET[as.character(snp)])
      class(res) <- "DNAbin"
      return(res)
    }
    
    ## create transversions for defined SNPs
    transv <- function(res){
      res <- sapply(TRANSVSET[as.character(res)],sample,1)
      class(res) <- "DNAbin"
      return(res)
    }
    
    ## duplicate a sequence (including possible mutations)
    seq.dupli <- function(res){
      n.transi <- rbinom(n=1, size=length.seqs, prob=mu.transi) # total number of transitions
      if(n.transi>0) {
        idx <- sample(length.seqs, size=n.transi, replace=FALSE)
        res[idx] <- transi(res[idx])
      }
      
      ## transversions ##
      n.transv <- rbinom(n=1, size=length.seqs, prob=mu.transv) # total number of transitions
      if(n.transv>0) {
        idx <- sample(length.seqs, size=n.transv, replace=FALSE)
        res[idx] <- transv(res[idx])
      }
      return(res)
    }
    
    class(res) <- "DNAbin"
    
    res <- replicate(num.seqs, seq.dupli(res))
    res <- as.matrix.DNAbin(res)
    
    class(res) <- "DNAbin"
    
    res <- t(res)
    
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
  
    if (!is.null(prop.haps) && prop.haps <= 1 / Hstar) {
      stop("prop.haps must be greater than 1 / Hstar")
    }
  
    if (!is.null(prop.seqs) && prop.seqs <= 1 / N)  {
      stop("prop.seqs must be greater than 1 / N")
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
	# individuals' haplotypes from the probabilities ##
	  
	  gen.perms <- function() {
	    if (is.null(subset.haps)) {
	      sample(haps, size = num.specs, replace = TRUE, prob = probs)
	    } else {
	      sample(subset.haps, size = num.specs, replace = TRUE, prob = probs[subset.haps])
	    }
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
	
	 if (is.null(subset.haps)) {
	   Q <- Hstar - P
	   assign("R", P / Hstar, envir = .GlobalEnv)
	   S <- (Hstar - P) / Hstar
	   assign("Nstar", (N * Hstar) / P, envir = .GlobalEnv)
	   assign("X", ((N * Hstar) / P) - N, envir = .GlobalEnv)
	 } else {
	   Q <- length(subset.haps) - P
	   assign("R", P / length(subset.haps), envir = .GlobalEnv)
	   S <- (length(subset.haps) - P) / length(subset.haps)
	   assign("Nstar", (N * length(subset.haps)) / P, envir = .GlobalEnv)
	   assign("X", ((N * length(subset.haps)) / P) - N, envir = .GlobalEnv)
	   }
	 
	## Calculate slope of curve using last n points (or proportion of points) on curve
	
	# perms must be large enough to ensure monotonicity and a non-negative slope
    
	  if (!is.null(prop.pts) && is.null(num.pts)) { 
	    lin.reg <- lm(means ~ specs, data = tail(d, n = ceiling(prop.pts * length(d))))
	    assign("beta1", abs(coef(lin.reg)[[2]]), envir = .GlobalEnv)
	  }
	 
	 if (!is.null(num.pts) && is.null(prop.pts)) {
	    lin.reg <- lm(means ~ specs, data = tail(d, n = num.pts))
	    assign("beta1", abs(coef(lin.reg)[[2]]), envir = .GlobalEnv)
	 }
	 
  ## Output results to R console and text file ##
		    
	  cat("\n \n --- Measures of Sampling Closeness --- \n \n", 
	  "Mean number of haplotypes sampled: " , P, 
	  "\n Mean number of haplotypes not sampled: " , Q, 
	  "\n Proportion of haplotypes (specimens) sampled: " , R, 
	  "\n Proportion of haplotypes (specimens) not sampled: " , S,
	  "\n \n Mean value of N*: ", Nstar / K,
	  "\n Mean number of specimens not sampled: ", X / K, 
	  "\n \n Haplotype accumulation curve slope: ", beta1,
	  "\n Mean number of specimens required to observe one new haplotype: ", 1 / beta1)
	  
	  cat("\n \n 95% CI for number of haplotypes recovered: ", c(max(lower), max(upper)))
	
    df[nrow(df) + 1, ] <- c(P, max(lower), max(upper), Q, R, S, Nstar / K, X / K, beta1, 1 / beta1)
    
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
	    abline(h = R * Hstar, v = N / K, lty = 2) # dashed line
	    abline(h = p * Hstar, lty = 3) # dotted line
	    HAC.bar <- barplot(num.specs * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = haps, main = "Haplotype frequency distribution")
	  } else {
	    abline(h = c(R * length(subset.haps)), v = N / K, lty = 2) # dashed line
	    abline(h = p * length(subset.haps), lty = 3) # dotted line
	    HAC.bar <- barplot(num.specs * (probs[subset.haps] / sum(probs[subset.haps])), xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = subset.haps, main = "Haplotype frequency distribution")
	  }
	  df
}
