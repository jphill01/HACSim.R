### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function() {
  iters <- 1
  HAC.sim(N = N, Hstar = Hstar, probs = probs, K = K, perms = perms, p = p, input.seqs = input.seqs)
	while (R < p) {
		HAC.sim(N = ceiling(Nstar), Hstar = Hstar, probs = probs, K = K, perms = perms, p = p)
	  iters <- iters + 1
	}
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("Desired level of H* has not yet been reached \n")
  } else{
    cat("Desired level of H* has been reached. \n \n \n The algorithm converged after", iters, "iterations.", "\n")
  }
}
