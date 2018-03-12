HAC.simest <- function(model = c("GAM", "SCAM", "Krig"), k = 10){
	
	if (model == "GAM") {
		
		cat("\n Thin plate smooth (tp) \n")
		HAC.tp <- inv.predict(HAC.tp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.tp)
				
		cat("\n Cubic spline smooth (cr) \n")
		HAC.cr <- inv.predict(HAC.cr, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.cr)
	
		cat("\n P-spline smooth (ps) \n")
		HAC.ps <- inv.predict(HAC.ps, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.ps)
	
		cat("\n Adaptive smooth (ad) \n")
		HAC.ad <- inv.predict(HAC.ad, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.ad)
	
	}
	
	if (model == "SCAM") {
		
		cat("\n Monotonically increasing smooth (mpi) \n")
		HAC.mpi <- inv.predict(HAC.mpi, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.mpi)
		
		cat("\n Concave smooth (cv) \n")
		HAC.cv <- inv.predict(HAC.cv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.cv)
	
		cat("\n Monotonically increasing and concave smooth (micv) \n")
		HAC.micv <- inv.predict(HAC.micv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.micv)
	
	}
	
	if (model == "Krig") {
		
		cat("\n Matern covariance function (Matern) \n")
		HAC.matern <- inv.predict(HAC.matern, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.matern)
		
		cat("\n Spherical covariance function (sph) \n")
		HAC.sph <- inv.predict(HAC.sph, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.sph)	
		
		cat("\n Exponential covariance function (exp) \n")
		HAC.exp <- inv.predict(HAC.exp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		print(HAC.exp)
		
	}
	
}
