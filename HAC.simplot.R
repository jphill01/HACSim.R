HAC.simplot <- function(model = c("GAM", "SCAM", "Krig"), k = 10) {
	
	if (model == "GAM") {
		
		cat("\n Thin plate smooth (tp) \n")
		HAC.tp <- gam(means ~ s(specs, bs = "tp", k = k), optimizer = c("outer", "bfgs"), data = d)
		gam.check(HAC.tp)
	
			cat("\n P-spline smooth (ps) \n")
		HAC.ps <- gam(means ~ s(specs, bs = "ps", k = k), optimizer = c("outer", "bfgs"), data = d)
		gam.check(HAC.ps)
		
		cat("\n Cubic spline smooth (cr) \n")
		HAC.cr <- gam(means ~ s(specs, bs = "cr", k = k), optimizer = c("outer", "bfgs"), data = d)
		gam.check(HAC.cr)
		
		cat("\n Adaptive smooth (ad) \n")
		HAC.ad <- gam(means ~ s(specs, bs = "ad", k = k), optimizer = c("outer", "bfgs"), data = d)
		gam.check(HAC.ad)

	} 
	
	
		
	if (model == "SCAM") {
		
		cat("\n Monotonically increasing smooth (mpi) \n")
		HAC.mpi <- scam(means ~ s(specs, bs = "mpi", k = k), data = d)
		scam.check(HAC.mpi)	
		
		cat("\n Concave smooth (cv) \n")
		HAC.cv <- scam(means ~ s(specs, bs = "cv", k = k), data = d)
		scam.check(HAC.cv)	
		
		cat("\n Monotonically increasing and concave smooth (micv) \n")
		HAC.micv <- scam(means ~ s(specs, bs = "micv", k = k), data = d)
		scam.check(HAC.micv)	
	}
	
	if (model == "Krig") {
	  
	  cat("\n Matern covariance function (Matern) \n")
	  HAC.matern <- gam(means ~ s(specs, bs = "gp", k = k), optimizer = c("outer", "bfgs"), data = d)
	  gam.check(HAC.matern)
		
		cat("\n Spherical covariance function (sph) \n")
		HAC.sph <- gam(means ~ s(specs, bs = "gp", k = k, m = 1), optimizer = c("outer", "bfgs"), data = d)
		gam.check(HAC.sph)
		
		cat("\n Exponential covariance function (exp) \n")
		HAC.exp <- gam(means ~ s(specs, bs = "gp", k = k, m = 2), optimizer = c("outer", "bfgs"), data = d)
		gam.check(HAC.exp)

	}

}