HAC.simest <- function(model = c("GAM", "SCAM", "Krig"), k = 10, method = c("Bisect", "Newton", "Both")){
	
	if (model == "GAM") {
		
		cat("\n Thin plate smooth (tp) \n")
	  if (method == "Bisect") {
	    HAC.tp.bisect <- inv.predict(HAC.tp, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
	    cat(HAC.tp.bisect, "\n")
	  }
	  if (method == "Newton") {
	    HAC.tp.newton <- inv.predict.newton(HAC.tp, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
	    cat(HAC.tp.newton, "\n")
	  }
    if (method == "Both") {
      HAC.tp.bisect <- inv.predict(HAC.tp, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
      cat(HAC.tp.bisect, "(Bisection) \n")
      HAC.tp.newton <- inv.predict.newton(HAC.tp, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.tp.newton, "(Newton) \n")
    }
				
		cat("\n Cubic spline smooth (cr) \n")
		if (method == "Bisect") {
		  HAC.cr.bisect <- inv.predict(HAC.cr, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.cr.bisect, "\n")
		}
		if (method == "Newton") {
		  HAC.cr.newton <- inv.predict.newton(HAC.cr, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.cr.newton, "\n")
		}
		if (method == "Both") {
		  HAC.cr.bisect <- inv.predict(HAC.cr, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.cr.bisect, "(Bisection) \n")
		  HAC.cr.newton <- inv.predict.newton(HAC.cr, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.cr.newton, "(Newton) \n")
		}
	
		cat("\n P-spline smooth (ps) \n")
		if (method == "Bisect") {
		  HAC.ps.bisect <- inv.predict(HAC.ps, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.ps.bisect, "\n")
		}
    if (method == "Newton") {
      HAC.ps.newton <- inv.predict.newton(HAC.ps, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.ps.newton, "\n")
    }
    if(method == "Both") {
      HAC.ps.bisect <- inv.predict(HAC.ps, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
      cat(HAC.ps.bisect, "(Bisection) \n")
      HAC.ps.newton <- inv.predict.newton(HAC.ps, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.ps.newton, "(Newton) \n")
    }
	
		cat("\n Adaptive smooth (ad) \n")
		if (method == "Bisect") {
		  HAC.ad.bisect <- inv.predict(HAC.ad, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.ad.bisect, "\n")
		}
		if (method == "Newton") {
		  HAC.ad.newton <- inv.predict.newton(HAC.ad, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.ps.newton, "\n")
		}
    if (method == "Both") {
      HAC.ad.bisect <- inv.predict(HAC.ad, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
      cat(HAC.ad.bisect, "(Bisection) \n")
      HAC.ad.newton <- inv.predict.newton(HAC.ad, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.ps.newton, "(Newton) \n")
    }
	
	}
	
	if (model == "SCAM") { 
	  
	  # bisection will work fine for SCAMs since the fit is forced to be monotonic in the region of interest
	  # the bisection method is problematic if fits are not monotonic or prediction is needed near the asymptote
		
		cat("\n Monotonically increasing smooth (mpi) \n")
	  if (method == "Bisect") {
	    HAC.mpi.bisect <- inv.predict(HAC.mpi, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
	    cat(HAC.mpi.bisect, "\n")
	  }
		if (method == "Newton") {
		  HAC.mpi.newton <- inv.predict.newton(HAC.mpi, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.mpi.newton, "\n")
		}
    if (method == "Both") {
      HAC.mpi.bisect <- inv.predict(HAC.mpi, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
      cat(HAC.mpi.bisect, "(Bisection) \n")
      HAC.mpi.newton <- inv.predict.newton(HAC.mpi, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.mpi.newton, "(Newton) \n")
    }
		
		cat("\n Concave smooth (cv) \n")
		if (method == "Bisect") {
		  HAC.cv.bisect <- inv.predict(HAC.cv, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.cv.bisect, "\n")
		}
		if (method == "Newton") {
		  HAC.cv.newton <- inv.predict.newton(HAC.cv, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.cv.newton, "\n")
		}
	  if (method == "Both") {
	    HAC.cv.bisect <- inv.predict(HAC.cv, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
	    cat(HAC.cv.bisect, "(Bisection) \n")
	    HAC.cv.newton <- inv.predict.newton(HAC.cv, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
	    cat(HAC.cv.newton, "(Bisection) \n")
	  }
		
		cat("\n Monotonically increasing and concave smooth (micv) \n")
		if (method == "Bisect") {
		  HAC.micv.bisect <- inv.predict(HAC.micv, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.micv.bisect, "\n")
		}
		if (method == "Newton") {
		  HAC.micv.newton <- inv.predict.newton(HAC.micv, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.micv.newton, "\n")
		}
    if (method == "Both") {
      HAC.micv.bisect <- inv.predict(HAC.micv, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
      cat(HAC.micv.bisect, "(Bisection) \n")
      HAC.micv.newton <- inv.predict.newton(HAC.micv, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.micv.newton, "(Newton) \n")
    }
	
	}
	
	if (model == "Krig") {
		cat("\n Matern covariance function (Matern) \n")
	  if (method == "Bisect") {
	    HAC.matern.bisect <- inv.predict(HAC.matern, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
	    cat(HAC.matern.bisect, "\n")
	  }
	  if (method == "Newton") {
	    HAC.matern.newton <- inv.predict.newton(HAC.matern, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
	    cat(HAC.matern.newton, "\n")
	  }
		if (method == "Both") {
		  HAC.matern.bisect <- inv.predict(HAC.matern, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.matern.bisect, "(Bisection) \n")
		  HAC.matern.newton <- inv.predict.newton(HAC.matern, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.matern.newton, "(Newton) \n")
		}
		
		cat("\n Spherical covariance function (sph) \n")
		if (method == "Bisect") {
		  HAC.sph.bisect <- inv.predict(HAC.sph, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.sph.bisect, "\n")
		}
		if (method == "Newton"){
		  HAC.sph.newton <- inv.predict.newton(HAC.sph, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.sph.newton, "\n")
		}
		if (method == "Both") {
		  HAC.sph.bisect <- inv.predict(HAC.sph, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.sph.bisect, "{Bisection) \n")
		  HAC.sph.newton <- inv.predict.newton(HAC.sph, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.sph.newton, "(Newton) \n")
		}
		
		cat("\n Exponential covariance function (exp) \n")
		if (method == "Bisect") {
		  HAC.exp.bisect <- inv.predict(HAC.exp, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.exp.bisect, "\n")
		}
    if (method == "Newton") {
      HAC.exp.newton <- inv.predict.newton(HAC.exp, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
      cat(HAC.exp.newton, "\n")
    }
		if (method == "Both") {
		  HAC.exp.bisect <- inv.predict(HAC.exp, y = R*Hstar, x.name = "specs", lower = 1, upper = max(d$specs), interval = TRUE)
		  cat(HAC.exp.bisect, "(Bisection) \n")
		  HAC.exp.newton <- inv.predict.newton(HAC.exp, y = R*Hstar, x.name = "specs", start = N, interval = TRUE)
		  cat(HAC.exp.newton, "(Newton) \n")
		}
	}
	
}
