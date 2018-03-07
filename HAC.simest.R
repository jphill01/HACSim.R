HAC.simest <- function(model = c("GAM", "SCAM", "Krig"), k = 10){
	
	if (model == "GAM") {
		
		cat("\n Thin plate smooth (tp) \n")
		HAC.tp <- gam(means ~ s(specs, bs = "tp", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.tp <- inv.predict(HAC.tp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.tp)
				
		cat("\n Cubic spline smooth (cr) \n")
		HAC.cr <- gam(means ~ s(specs, bs = "cr", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.cr <- inv.predict(HAC.cr, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.cr)
	
		cat("\n P-spline smooth (ps) \n")
		HAC.ps <- gam(means ~ s(specs, bs = "ps", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.ps <- inv.predict(HAC.ps, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.ps)
	
			cat("\n Adaptive smooth (ad) \n")
		HAC.ad <- gam(means ~ s(specs, bs = "ad", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.ad <- inv.predict(HAC.ad, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.ad)
	
	}
	
	if (model == "SCAM") {
		
			cat("\n Monotonically increasing smooth (mpi) \n")
		HAC.mpi <- scam(means ~ s(specs, bs = "mpi", k = k), data = d)
		HAC.mpi <- inv.predict(HAC.mpi, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.mpi)
		
			cat("\n Concave smooth (cv) \n")
		HAC.cv <- scam(means ~ s(specs, bs = "cv", k = k), data = d)
		HAC.cv <- inv.predict(HAC.cv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.cv)
	
		cat("\n Monotonically increasing and concave smooth (micv) \n")
		HAC.micv <- scam(means ~ s(specs, bs = "micv", k = k), data = d)
		HAC.micv <- inv.predict(HAC.micv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.micv)
	
	}
	
	if (model == "Krig") {
		
		cat("\n Matern covariance function (Matern) \n")
		HAC.matern <- gam(means ~ s(specs, bs = "gp", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.matern <- inv.predict(HAC.matern, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.matern)
		
		cat("\n Spherical covariance function (sph) \n")
		HAC.sph <- gam(means ~ s(specs, bs = "gp", k = k, m = 1), optimizer = c("outer", "bfgs"), data = d)
		HAC.sph <- inv.predict(HAC.sph, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.sph)	
		
		cat("\n Exponential covariance function (exp) \n")
		HAC.exp <- gam(means ~ s(specs, bs = "gp", k = k, m = 2), optimizer = c("outer", "bfgs"), data = d)
		HAC.exp <- inv.predict(HAC.exp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.exp)
		
	}
	
}
