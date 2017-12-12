# HAC Simulation Iteration

HAC.simrep <- function(){
	while (R < p){
		HAC.sim(K = K, N = ceiling(Nstar), Hstar = Hstar, probs = probs, m = m, perms = perms, p = p)
	}
}