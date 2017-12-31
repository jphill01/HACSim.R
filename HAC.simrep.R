### HAC Simulation Iteration ###

## Run HAC Simulator until convergence is reached ##

HAC.simrep <- function() {
	while (R < p) {
		HAC.sim(N = ceiling(Nstar), Hstar = Hstar, probs = probs, K = K, perms = perms, p = p, plot.out = plot.out)
	}
}



