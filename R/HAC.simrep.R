### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function() {
  ptm <- proc.time()
  iters <- 1
    HAC.sim(N = N, 
            Hstar = Hstar, 
            probs = probs, 
            perms = perms, 
            p = p,
            input.seqs = input.seqs,
            subset.seqs = subset.seqs,
            subset.haps = subset.haps,
            prop.seqs = prop.seqs,
            prop.haps = prop.haps
            )
	  while (R < p) {
		  HAC.sim(N = ceiling(Nstar), 
		          Hstar = Hstar, 
		          probs = probs,
		          perms = perms,
		          p = p,
		          subset.seqs = subset.seqs,
		          subset.haps = subset.haps,
		          prop.seqs = prop.seqs,
		          prop.haps = prop.haps
            )
	    iters <- iters + 1
	    amt <- proc.time() - ptm
    }
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("Desired level of H* has not yet been reached \n")
  } else {
      cat("\n Desired level of H* has been reached. \n \n The algorithm converged after", 
        iters, "iterations and took", amt[3], "s.", "\n")
  }

}
