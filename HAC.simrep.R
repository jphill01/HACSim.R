### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

iters <- 1

HAC.simrep <- function() {
  HAC.sim(N = N, Hstar = Hstar, probs = probs, K = K, perms = perms, p = p, input.seqs = input.seqs)
	while (R < p) {
		assign("iters", iters + 1, envir = .GlobalEnv)
		HAC.sim(N = ceiling(Nstar), Hstar = Hstar, probs = probs, K = K, perms = perms, p = p)
	}
}
