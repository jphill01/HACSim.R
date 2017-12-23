### Haplotype Accumulation Curve Simulator ###

HAC.sim <- function(K = 1, N, Hstar, probs, perms = 10000, p = 1, plot.out = TRUE) {
	
	## Error messages ##
	
	if (N > H) {
		stop("N must be greater than or equal to H")
	}
	
	if (N > K) {
		stop("N must be greater than or equal to K")
	}
	
	if (sum(probs) != 1) {
		stop("probs must sum to 1")
	}
	
	## Set up container(s) to hold the identity of each individual from each permutation ##
	
	num.specs <- ceiling(N / K)
	
	pop <- array(dim = c(c(perms, num.specs), K))
	
	## Create an ID for each haplotype ##
	
	haps <- as.character(1:Hstar)
	
	## Assign individuals (N) to each subpopulation (K) ##
	
	specs <- 1:num.specs
	
	## Generate permutations, assume each permutation has N individuals, and sample those individuals' haplotypes from the probabilities ##
	
	for (j in 1:perms) {
		for (i in 1:K) {
				pop[j, specs, i] <- sample(haps, size = num.specs, replace = TRUE, prob = probs)
			}
	}
	
	## Make a matrix to hold individuals from each permutation ##

	HAC.mat <- array(dim = c(c(perms, num.specs), K))
	
	## Perform haplotype accumulation ##
	
	for (k in specs) {
		for (j in 1:perms) {
			for (i in 1:K) {
				ind.index <- sample(specs, size = k, replace = FALSE) # which individuals are sampled
				hap.plot <- pop[sample(1:nrow(pop), size = 1, replace = TRUE), ind.index, sample(i, size = 1, replace = TRUE)] # extract those individuals from a permutation
				HAC.mat[j, k, i] <- length(unique(hap.plot)) # how many haplotypes recovered a given sampling intensity (k) from each permutation (j)
			}
		}
	}
	
	## Calculate the mean and CI for number of haplotypes recovered

	means <- apply(HAC.mat, MARGIN = 2, mean)
	lower <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.025))
	upper <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.975))
	
	## Make data accessible to user ##

	assign("d", data.frame(specs, means), envir = .GlobalEnv)
	
	## Compute simple summary statistics and display output ##
	
	P <- max(means)
	Q <- Hstar - max(means)
	assign("R", max(means) / Hstar, envir = .GlobalEnv)
	S <- (Hstar - max(means)) / Hstar
	assign("Nstar", (N * Hstar) / max(means), envir = .GlobalEnv)
	X <- ((N * Hstar) / max(means)) - N
				    
	cat("\n Measures of Sampling Closeness \n \n Mean number of haplotypes sampled: " , P, "\n Mean number of haplotypes not sampled: " , Q, "\n Proportion of haplotypes sampled: " , R, "\n Proportion of haplotypes not sampled:  " , S, "\n \n Mean value of N*: ", Nstar, "\n Mean number of individuals not sampled: ", X, "\n \n")
	
	## Check whether desired level of haplotype recovery has been reached ##
		
	if (R < p) {
		cat("Desired level of H* has not yet been reached \n")
		} else{
			cat("Desired level of H* has been reached")
	}
	
	## Plot the haplotype accumulation curve and haplotype frequency barplot ##
	
	if (plot.out == TRUE) {
		par(mfrow = c(1, 2))
		plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes",  ylim = c(1, Hstar))
		polygon(x = c(specs, rev(specs)), y = c(lower, rev(upper)), col = "gray")
		lines(specs, means, lwd = 2)
		HAC.bar <- barplot(num.specs * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = 1:Hstar)
		}

}
