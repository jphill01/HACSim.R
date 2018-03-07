##### Best Model Output

summary.HAC.sim <- function() {
  
		#cat("\n Thin plate smooth (tp) \n")
		HAC.tp1 <- gam(means ~ s(specs, bs = "tp", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.tp2 <- inv.predict(HAC.tp1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.tp)
				
		#cat("\n Cubic spline smooth (cr) \n")
		HAC.cr1 <- gam(means ~ s(specs, bs = "cr", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.cr2 <- inv.predict(HAC.cr1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.cr)
	
		#cat("\n P-spline smooth (ps) \n")
		HAC.ps1 <- gam(means ~ s(specs, bs = "ps", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.ps2 <- inv.predict(HAC.ps1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.ps)
	
		#cat("\n Adaptive smooth (ad) \n")
		HAC.ad1 <- gam(means ~ s(specs, bs = "ad", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.ad2 <- inv.predict(HAC.ad1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.ad)
		

		#cat("\n Monotonically increasing smooth (mpi) \n")
		HAC.mpi1 <- scam(means ~ s(specs, bs = "mpi", k = k), data = d)
		HAC.mpi2 <- inv.predict(HAC.mpi1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.micv)
		
		#cat("\n Concave smooth (cv) \n")
		HAC.cv1 <- scam(means ~ s(specs, bs = "cv", k = k), data = d)
		HAC.cv2 <- inv.predict(HAC.cv1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.cv)
	
		#cat("\n Monotonically increasing and concave smooth (micv) \n")
		HAC.micv1 <- scam(means ~ s(specs, bs = "micv", k = k), data = d)
		HAC.micv2 <- inv.predict(HAC.micv1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.micv)
	

		#cat("\n Matern covariance function (Matern) \n")
		HAC.matern1 <- gam(means ~ s(specs, bs = "gp", k = k), optimizer = c("outer", "bfgs"), data = d)
		HAC.matern2 <- inv.predict(HAC.matern1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.matern)
		
		#cat("\n Spherical covariance function (sph) \n")
		HAC.sph1 <- gam(means ~ s(specs, bs = "gp", k = k, m = 1), optimizer = c("outer", "bfgs"), data = d)
		HAC.sph2 <- inv.predict(HAC.sph1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.sph)	
		
		#cat("\n Exponential covariance function (exp) \n")
		HAC.exp1 <- gam(means ~ s(specs, bs = "gp", k = k, m = 2), optimizer = c("outer", "bfgs"), data = d)
		HAC.exp2 <- inv.predict(HAC.exp1, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
		#print(HAC.exp)
		
			
	cat("----- Model Summary -----
    \n Best model: ",
		"\n AIC: ", min(c(HAC.tp1$aic, HAC.cr1$aic, HAC.ps1$aic, HAC.ad1$aic, HAC.mpi1$aic, HAC.cv1$aic, HAC.micv1$aic,  HAC.matern1$aic, HAC.sph1$aic, HAC.exp1$aic)),
		"\n Bootstrap SE: ",
		"\n Bootstrap bias: ",
		"\n Bootstrap % relative bias: ",
		"\n Bootstrap 95% CI: "
		)
}