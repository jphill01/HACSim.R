steff <- function() {
  N1 <- HAC.sim(N = ceiling(N), Hstar = Hstar, probs = probs, K = K, perms = perms, p = p)
  N2 <- HAC.sim(N = ceiling(Nstar), Hstar = Hstar, probs = probs, K = K, perms = perms, p = p)
		while (R < p) {
			N3 <- N - ((N1 - N)^2 / (N2 - 2 * N1 + N))
			N <- N3
			
	} 
}