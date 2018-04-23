HAC.simest <- function(model = c("GAM", "SCAM", "Krig"), k = 10){
	
	if (model == "GAM") {
		
		cat("\n Thin plate smooth (tp) \n")
    HAC.tp.bisect <- inv.predict(HAC.tp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.tp.bisect, "(bisection) \n")
		HAC.tp.newton <- inv.predict.newton(HAC.tp, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		cat(HAC.tp.newton, "(Newton) \n \n")
				
		cat("\n Cubic spline smooth (cr) \n")
    HAC.cr.bisect <- inv.predict(HAC.cr, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.cr.bisect, "(bisection) \n")
		HAC.cr.newton <- inv.predict.newton(HAC.cr, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		cat(HAC.cr.newton, "(Newton) \n \n")
	
		cat("\n P-spline smooth (ps) \n")
    HAC.ps.bisect <- inv.predict(HAC.ps, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.ps.bisect, "(bisection) \n")
		HAC.ps.newton <- inv.predict.newton(HAC.ps, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		cat(HAC.ps.newton, "(Newton) \n \n")
	
		cat("\n Adaptive smooth (ad) \n")
    HAC.ad.bisect <- inv.predict(HAC.ad, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.ad.bisect, "(bisection) \n")
		HAC.ad.newton <- inv.predict.newton(HAC.ad, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		cat(HAC.ps.newton, "(Newton) \n \n")
	
	}
	
	if (model == "SCAM") { 
	  
	  # bisection will work fine for SCAMs since the fit is forced to be monotonic in the region of interest
	  # the bisection method is problematic if fits are not monotonic or predction is needed near the asymptote
		
		cat("\n Monotonically increasing smooth (mpi) \n")
		HAC.mpi.bisect <- inv.predict(HAC.mpi, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
		cat(HAC.mpi.bisect, "(bisection) \n")
    HAC.mpi.newton <-inv.predict.newton(HAC.mpi, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
    cat(HAC.mpi.newton, "(Newton) \n \n")
		
		cat("\n Concave smooth (cv) \n")
		HAC.cv.bisect <- inv.predict(HAC.cv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
		cat(HAC.cv.bisect, "(bisection) \n")
    HAC.cv.newton <- inv.predict.newton(HAC.cv, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
    cat(HAC.cv.newton, "(Newton) \n \n")
	
		cat("\n Monotonically increasing and concave smooth (micv) \n")
		HAC.micv.bisect <- inv.predict(HAC.micv, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
		cat(HAC.micv.bisect, "(bisection) \n")
    HAC.micv.newton <- inv.predict.newton(HAC.micv, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
    cat(HAC.micv.newton, "(Newton) \n \n")
	
	}
	
	if (model == "Krig") {
		
		cat("\n Matern covariance function (Matern) \n")
    HAC.matern.bisect <- inv.predict(HAC.matern, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.matern.bisect, "(bisection) \n")
		HAC.matern.newton <- inv.predict.newton(HAC.matern, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
    cat(HAC.matern.newton, "(Newton) \n \n")
		
		cat("\n Spherical covariance function (sph) \n")
    HAC.sph.bisect <- inv.predict(HAC.sph, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.sph.bisect, "(bisection) \n")
		HAC.sph.newton <- inv.predict.newton(HAC.sph, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		cat(HAC.sph.newton, "(Newton) \n \n")
		
		cat("\n Exponential covariance function (exp) \n")
    HAC.exp.bisect <- inv.predict(HAC.exp, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = TRUE)
    cat(HAC.exp.bisect, "(bisection) \n")
		HAC.exp.newton <- inv.predict.newton(HAC.exp, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		cat(HAC.exp.newton, "(Newton) \n \n")
	}
	
}
