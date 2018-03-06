### Semiparametric Model Fitting ###

HAC.simfit <- function(model = c("GAM", "SCAM", "Krig"), k = 10) {
	
	if (model == "GAM") {
		
		cat("\n \n Thin plate smooth (tp) \n")
		HAC.tp <- gam(means ~ s(specs, bs = "tp", k = k), optimizer = c("outer", "bfgs"), data = d)
		cat("\n AIC: " ,  HAC.tp$aic, "\n")
		
		cat("\n \n Cubic spline smooth (cr) \n")
		HAC.cr <- gam(means ~ s(specs, bs = "cr", k = k), optimizer = c("outer", "bfgs"), data = d)
		cat("\n AIC: " ,  HAC.cr$aic, "\n")
		
		cat("\n \n P-spline smooth (ps) \n")
		HAC.ps <- gam(means ~ s(specs, bs = "ps", k = k), optimizer = c("outer", "bfgs"), data = d)
		cat("\n AIC: " ,  HAC.ps$aic, "\n")
		
		cat("\n \n Adaptive smooth (ad) \n")
		HAC.ad <- gam(means ~ s(specs, bs = "ad", k = k), optimizer = c("outer", "bfgs"), data = d)
		cat("\n AIC: " ,  HAC.ad$aic, "\n")
	} 
		
	if (model == "SCAM") {
		
		cat("\n \n Monotonically increasing smooth (mpi) \n")	
		HAC.mpi <- scam(means ~ s(specs, bs = "mpi", k = k), data = d)
		cat("\n AIC: " ,  HAC.mpi$aic, "\n")
		
		cat("\n \n Concave smooth (cv) \n")	
		HAC.cv <- scam(means ~ s(specs, bs = "cv", k = k), data = d)
		cat("\n AIC: " ,  HAC.cv$aic, "\n")
		
		cat("\n \n Monotonically increasing and concave smooth (micv) \n")	
		HAC.micv <- scam(means ~ s(specs, bs = "micv", k = k), data = d)
		cat("\n AIC: " ,  HAC.micv$aic, "\n")
	}
	
	if (model == "Krig") {
	  
	  cat("\n Matern covariance function \n")
	  HAC.matern <- gam(means ~ s(specs, bs = "gp", k = k), optimizer = c("outer", "bfgs"), data = d)
	  cat("\n AIC: " ,  HAC.matern$aic, "\n")
		
		cat("\n Spherical covariance function (sph) \n")
		HAC.sph <- gam(means ~ s(specs, bs = "gp", k = k, m = 1), optimizer = c("outer", "bfgs"), data = d)
		cat("\n AIC: " ,  HAC.sph$aic, "\n")
		
		cat("\n Exponential covariance function \n")
		HAC.exp <- gam(means ~ s(specs, bs = "gp", k = k, m = 2), optimizer = c("outer", "bfgs"), data = d)
		cat("\n AIC: " ,  HAC.exp$aic, "\n")
		
	}

}
