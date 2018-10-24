HAC.simfit <- function(model = c("GAM", "SCAM", "Krig"), k = 10) {
	
	if (model == "GAM") {
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.tp, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Thin plate smooth")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.cr, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Cubic smooth")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.ps, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "P-spline smooth")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.ad, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Adaptive smooth")
		lines(xx, yy, lty = 2, col = "red")
		
	}
	
		if (model == "SCAM") {
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.mpi, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Monotonically increasing smooth")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.cv, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Concave smooth")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.micv, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Monotonically increasing and concave smooth")
		lines(xx, yy, lty = 2, col = "red")
		
	}	

		if (model == "Krig") {
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.matern, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Matern covariance function")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.sph, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Spherical covariance function")
		lines(xx, yy, lty = 2, col = "red")
		
		xx <- seq(from = min(d$specs), to = max(d$specs))
		yy <- predict(HAC.exp, newdata = data.frame(specs = xx))
		plot(d$specs, d$means, main = "Exponential covariance function")
		lines(xx, yy, lty = 2, col = "red")
		
	}	
	
}

