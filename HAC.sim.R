### Haplotype Accumulation Curve Simulator ###

HAC.sim <- function(N, Hstar, probs, K = 1, m = 0, perms = 10000, p = 1, plot.out = TRUE) {
	
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
	
	if (sum(probs) != 1){
		stop("probs must sum to 1")
	}
	
	## Set up container(s) to hold the identity of each individual from each permutation ##

	num.specs <- ceiling(N / K)
	
	## Create an ID for each haplotype ##
	
	haps <- 1:Hstar
	
	## Assign individuals (N) to each subpopulation (K) ##
	
	specs <- 1:num.specs
	
	## Generate permutations, assume each permutation has N individuals, and sample those individuals' haplotypes from the probabilities ##
	
	gen.perms <- function() {
		sample(haps, size = num.specs, replace = TRUE, prob = probs)
	}
	
	pop <- array(dim = c(perms, num.specs, K))
	
	for (i in 1:K) {
		pop[,, i] <- replicate(perms, gen.perms())
	}
	
	## Allow individuals to migrate between subpopulations according to migration rate m ##
	
	if (m != 0) {
		for (i in 1:K) {
			for (j in 1:K) {
				i <- sample(perms, size = ceiling(perms * m), replace = FALSE)
				j <- sample(perms, size = ceiling(perms * m), replace = FALSE)
				tmp <- pop[i,, ]
				pop[i,, ] <- pop[j,, ]
				pop[j,, ] <- tmp
			}
		}
	}

	## Perform haplotype accumulation ##
	
	HAC.mat <- fillCube(pop, specs, perms, K)
	
	## Calculate the mean and CI for number of haplotypes recovered

	means <- apply(HAC.mat, MARGIN = 2, mean)
	lower <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.025))
	upper <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.975))
	
	## Make data accessible to user ##

	assign("d", data.frame(specs, means), envir = .GlobalEnv)
	
	## Compute simple summary statistics and display output ##
	
	P <- tail(means, n = 1)
	Q <- Hstar - tail(means, n = 1)
	assign("R", tail(means, n = 1) / Hstar, envir = .GlobalEnv)
	S <- (Hstar - tail(means, n = 1)) / Hstar
	assign("Nstar", (N * Hstar) / tail(means, n = 1), envir = .GlobalEnv)
	X <- ((N * Hstar) / tail(means, n = 1)) - N
				    
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
