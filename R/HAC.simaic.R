### Semiparametric Model Fitting ###

HAC.simaic <- function(model = c("GAM", "SCAM", "Krig"), k = 10) {
	
	if (model == "GAM") {
		
		cat("\n \n Thin plate smooth (tp) \n")
		cat("\n AIC: " ,  HAC.tp$aic, "\n")
		
		cat("\n \n Cubic spline smooth (cr) \n")
		cat("\n AIC: " ,  HAC.cr$aic, "\n")
		
		cat("\n \n P-spline smooth (ps) \n")
		cat("\n AIC: " ,  HAC.ps$aic, "\n")
		
		cat("\n \n Adaptive smooth (ad) \n")
		cat("\n AIC: " ,  HAC.ad$aic, "\n")
		
		}
		
		if (model == "SCAM") {
		
		cat("\n Monotonically increasing smooth (mpi) \n")	
		cat("\n AIC: " ,  HAC.mpi$aic, "\n")
		
		cat("\n \n Concave smooth (cv) \n")	
		cat("\n AIC: " ,  HAC.cv$aic, "\n")
		
		cat("\n \n Monotonically increasing and concave smooth (micv) \n")	
		cat("\n AIC: " ,  HAC.micv$aic, "\n")
		
		}
		
		if (model == "Krig") {
	  
	  cat("\n Matern covariance function (matern) \n")
	  cat("\n AIC: " ,  HAC.matern$aic, "\n")
		
		cat("\n Spherical covariance function (sph) \n")
		cat("\n AIC: " ,  HAC.sph$aic, "\n")
		
		cat("\n Exponential covariance function (exp) \n")
		cat("\n AIC: " ,  HAC.exp$aic, "\n")
		
		}

}
