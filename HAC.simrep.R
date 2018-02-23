### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function() {
	while (R < p) {
		HAC.sim(N = ceiling(Nstar), Hstar = Hstar, probs = probs, K = K, m = m, perms = perms, p = p, seqs = seqs)
	}
}



