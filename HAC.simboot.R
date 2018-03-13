### Bootstrap Simulation ###

HAC.simboot <- function(model = c("GAM", "SCAM", "Krig"), k = 10) {
  
  ## Set progress bar ##
  
  #if (progress == TRUE) {
    #pb <- utils::txtProgressBar(min = 0, max = iters, style = 3)
  #}
	
	if (model == "GAM"){
		
		cat("\n Thin plate smooth (tp) \n")
		res <- resid(HAC.tp) - mean(resid(HAC.tp)) # centre the residuals
		n <- length(res)
		boot.data <- data.frame(d, res = res, fit = fitted(HAC.tp))
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "tp", k = k), optimizer = c("outer", "bfgs"), data = data)
		  # Simulate the correct variance
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		# Make sure the original estimate also gets returned
		if (all(i == 1:n)) {
				inv.predict(HAC.tp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0 # relative bias
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))
		
		
		cat("\n Cubic spline smooth (cr) \n")
		res <- resid(HAC.cr) - mean(resid(HAC.cr)) 
		n <- length(res)
		boot.data <- data.frame(d, res = res, fit = fitted(HAC.cr))
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "cr", k = k), optimizer = c("outer", "bfgs"), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		if (all(i == 1:n)) {
			inv.predict(HAC.cr, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0 # relative bias
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))     
    
    
    cat("\n P-spline smooth (ps) \n")
		res <- resid(HAC.ps) - mean(resid(HAC.ps)) 
		n <- length(res)
		boot.data <- data.frame(d, res = res, fit = fitted(HAC.ps))
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "ps", k = k), optimizer = c("outer", "bfgs"), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		if (all(i == 1:n)) {
			inv.predict(HAC.ps, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))
		
		
		cat("\n Adaptive smooth (ad) \n")
		res <- resid(HAC.ad) - mean(resid(HAC.ad))
		n <- length(res)
		boot.data <- data.frame(d, res = res, fit = fitted(HAC.ad))
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "ad", k = k), optimizer = c("outer", "bfgs"), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		if (all(i == 1:n)) {
			inv.predict(HAC.ad, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))    

	}
    
    if (model == "SCAM"){
    	
    cat("\n Monotonically increasing smooth (mpi) \n")
		res <- resid(HAC.mpi) - mean(resid(HAC.mpi))  
		n <- length(res)
		boot.data <- data.frame(d,  fit = fitted(HAC.mpi), res = res)
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- scam(boot.data$means + res[i] ~ s(specs, bs = "mpi", k = k), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		if (all(i == 1:n)) {
			inv.predict(HAC.mpi, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all")) 
		
		
		cat("\n Concave smooth (cv) \n")
		res <- resid(HAC.cv) - mean(resid(HAC.cv))  
		n <- length(res)
		boot.data <- data.frame(d,  fit = fitted(HAC.cv), res = res)
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- scam(boot.data$means + res[i] ~ s(specs, bs = "cv", k = k), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		if (all(i == 1:n)) {
			inv.predict(HAC.cv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))  
 

		cat("\n Monotonically increasing and concave smooth (micv) \n")
		res <- resid(HAC.micv) - mean(resid(HAC.micv))  
		n <- length(res)
		boot.data <- data.frame(d,  fit = fitted(HAC.micv), res = res)
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  boot.fit <- scam(boot.data$means + res[i] ~ s(specs, bs = "micv", k = k), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
		if (all(i == 1:n)) {
			inv.predict(HAC.micv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}

		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100) 
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))  
 
	}
 
 	if (model == "Krig") {
 	  
 	  cat("\n Matern covariance function \n")
 	  res <- resid(HAC.matern) - mean(resid(HAC.matern))  
 	  n <- length(res)
 	  boot.data <- data.frame(d, res = res, fit = fitted(HAC.matern))
 	  boot.fun <- function(data, i) {
 	    
 	    ## Update progress bar ##
 	    
 	    #if (progress == TRUE) {
 	      #utils::setTxtProgressBar(pb, i)
 	    #}
 	    
 	    boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "gp", k = k), optimizer = c("outer", "bfgs"), data = data)
 	    Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
 	  if (all(i == 1:n)) {
 	    inv.predict(HAC.matern, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
 	    } else {
 	      inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
 	    }
 	  }
 	  
 	  res <- boot(boot.data, boot.fun, R = 1000)  
 	  print(res)
 	  rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
 	  cat("\n Percent relative bias: ", rel.bias*100)
 	  plot(res)
 	  cat("\n\n")
 	  print(boot.ci(res, type = "all"))  
 	  
 		
 		cat("\n Spherical covariance function (sph) \n")
 		res <- resid(HAC.sph) - mean(resid(HAC.sph))  
		n <- length(res)
		boot.data <- data.frame(d, res = res, fit = fitted(HAC.sph))
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
		  
 		  boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "gp", k = k, m = 1), optimizer = c("outer", "bfgs"), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
 		if (all(i == 1:n)) {
			inv.predict(HAC.sph, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}
 
 		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100)
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))  

	  cat("\n Exponential covariance function (exp) \n")
 		res <- resid(HAC.exp) - mean(resid(HAC.exp))  
		n <- length(res)
		boot.data <- data.frame(d, res = res, fit = fitted(HAC.exp))
		boot.fun <- function(data, i) {
		  
		  ## Update progress bar ##
		  
		  #if (progress == TRUE) {
		    #utils::setTxtProgressBar(pb, i)
		  #}
		  
 		  boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = "gp", k = k, m = 2), optimizer = c("outer", "bfgs"), data = data)
		  Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
 		if (all(i == 1:n)) {
				inv.predict(HAC.exp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
			} else {
				inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		}
}
 
 		res <- boot(boot.data, boot.fun, R = 1000)  
		print(res)
		rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0
		cat("\n Percent relative bias: ", rel.bias*100)
		plot(res)
		cat("\n\n")
		print(boot.ci(res, type = "all"))

	}

}