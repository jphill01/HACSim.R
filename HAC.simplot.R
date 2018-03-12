HAC.simplot <- function(model = c("GAM", "SCAM", "Krig"), k = 10) {
	
	if (model == "GAM") {
		
		cat("\n Thin plate smooth (tp) \n")
		gam.check(HAC.tp)
	
		cat("\n P-spline smooth (ps) \n")
		gam.check(HAC.ps)
		
		cat("\n Cubic spline smooth (cr) \n")
		gam.check(HAC.cr)
		
		cat("\n Adaptive smooth (ad) \n")
		gam.check(HAC.ad)

	} 
	
	
		
	if (model == "SCAM") {
		
		cat("\n Monotonically increasing smooth (mpi) \n")
		scam.check(HAC.mpi)	
		
		cat("\n Concave smooth (cv) \n")
		scam.check(HAC.cv)	
		
		cat("\n Monotonically increasing and concave smooth (micv) \n")
		scam.check(HAC.micv)	
	}
	
	if (model == "Krig") {
	  
	  cat("\n Matern covariance function (Matern) \n")
	  gam.check(HAC.matern)
		
		cat("\n Spherical covariance function (sph) \n")
		gam.check(HAC.sph)
		
		cat("\n Exponential covariance function (exp) \n")
		gam.check(HAC.exp)

	}

}