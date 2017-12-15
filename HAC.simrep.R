### HAC Simulation Iteration ###

## Run HAC Simulator until convergence is reached ##

HAC.simrep <- function(){
	while (R < p){
		HAC.sim(K = K, N = ceiling(Nstar), Hstar = Hstar, probs = probs, m = m, perms = perms, p = p, plot.out = = plot.out)
	}
}



