### Haplotype Accumulation Curve Simulator ###

HAC.sim <- function(K = 1, N, Hstar, probs, m = 0, perms = 10000, p = 0.95){

	pop <- array(dim = c(c(perms, ceiling((1 - m) * N / K)), K))
	haps <- as.character(1:Hstar)
	specs <- 1:ceiling((1 - m) * N / K)
	
	for (j in 1:perms){
		for (i in 1:K){
			pop[j, specs, i] <- sample(haps, size = length(specs), replace = TRUE, prob = probs)
		}
	}
	

	HAC.mat <- array(dim = c(c(perms, length(specs), K)))
	
	for (k in specs){
		for (j in 1:perms){
			for (i in 1:K){
				ind.index <- sample(specs, size = k, replace = FALSE) 
				hap.plot <- pop[sample(1:nrow(pop), size = 1, replace = TRUE), ind.index, sample(i, size = 1, replace = TRUE)] 
				HAC.mat[j, k, i] <- length(unique(hap.plot))
			}
		}
	}

	means <- apply(HAC.mat, MARGIN = 2, mean)
	lower <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.025))
	upper <- apply(HAC.mat, MARGIN = 2, function(x) quantile(x, 0.975))

	d <- assign("d", data.frame(specs, means), envir = .GlobalEnv)
	
    P <- max(means)
	Q <- Hstar - max(means)
	R <- assign("R", max(means) / Hstar, envir = .GlobalEnv)
	S <- (Hstar - max(means)) / Hstar
	Nstar <- assign("Nstar", (N * Hstar) / max(means), envir = .GlobalEnv)
	X <- ((N * Hstar) / max(means)) - N
		    
	cat("\n Measures of Sampling Closeness \n \n Mean number of haplotypes sampled: " , P, "\n Mean number of haplotypes not sampled: " , Q, "\n Proportion of haplotypes sampled: " , R, "\n Proportion of haplotypes not sampled:  " , S, "\n \n Calculated mean value of N*: ", Nstar, "\n Mean number of individuals not sampled: ", X, "\n \n")
		
	if (R < p){
		cat("Desired level of H* has not yet been reached \n")
		} else{
			cat("Desired level of H* has been reached")
	}

	par(mfrow = c(1, 2))

	plot(specs, means, type = "n", xlab = "Specimens sampled", ylab = "Unique haplotypes",  ylim = c(1, Hstar))
	polygon(x = c(specs, rev(specs)), y = c(lower, rev(upper)), col = "gray")
	lines(specs, means, lwd = 2)
	HAC.bar <- barplot(length(specs) * probs, xlab = "Unique haplotypes", ylab = "Specimens sampled", names.arg = 1:Hstar)

}
