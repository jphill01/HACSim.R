# SCAM Variance-Covariance Matrix

##' @export
`vcov.scam` <- function(object, freq = FALSE, dispersion = NULL,
                         parametrized = TRUE, ...)  {
  if (freq) {
    vc <- if (parametrized) {
      object$Ve.t
    } else {
      object$Ve
    }
  } else {
    vc <- if (parametrized) {
      object$Vp.t
    } else {
      object$Vp
    }
  }
  if (!is.null(dispersion)) {
    vc <- dispersion * vc/object$sig2
  }
  name <- names(object$edf)
  dimnames(vc) <- list(name, name)
  vc
}

##' Extract coefficients from a fitted `scam` model.
##'
##' @param object a model object fitted by `scam()`
##' @param parametrized logical; extract parametrized coefficients, which respect the linear inequality constraints of the model.
##' @param ... other arguments.
##'
##' @export
`coef.scam` <- function(object, parametrized = TRUE, ...) {
  coefs <- if (parametrized) {
    object$coefficients.t
  } else {
    object$coefficients
  }
  coefs
}

# MVN Distribution Posterior Simulation

HAC.simpost <- function(model = c("GAM", "SCAM", "Krig"), k = 10){
	
	n <- 1000
	g <- data.frame(specs = seq(1, Nstar - X, length = 1000)) # fine grid for predictor

	if (model == "GAM") {
	  
	  cat("\n Thin plate smooth \n\n")
	  gg  <- predict(HAC.tp, newdata = g, type = "response")
	  Xp <- predict(HAC.tp, g, type = "lpmatrix") 
	  beta <- coef(HAC.tp) 
	  Vb <- vcov(HAC.tp)
	  mrand <- mvrnorm(n, beta, Vb)
	  opt <- rep(NA, n)
	  ilink <- family(HAC.tp)$linkinv
	  for (i in seq_len(n)) {
	    pred <- ilink(Xp %*% mrand[i, ])
	    opt[i] <- g$specs[which.max(pred)]
	  }
	  
      mu <- mean(opt) # mean
      cat("N*: ", mu, "\n")
      
      se <- sd(opt) / sqrt(n) # standard error
      cat("SE: ", se , "\n")
      
      ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
      cat("95% CI: ", ci , "\n")
      
      cat(opt)
      cat(quantile(opt, c(0.025, 0.975)))
	  
	  cat("\n Cubic spline smooth \n\n")
	  gg  <- predict(HAC.cr, newdata = g, type = "response")
	  Xp <- predict(HAC.cr, g, type = "lpmatrix") 
	  beta <- coef(HAC.cr) 
	  Vb <- vcov(HAC.cr)
	  mrand <- mvrnorm(n, beta, Vb)
	  opt <- rep(NA, n)
	  ilink <- family(HAC.cr)$linkinv
	  for (i in seq_len(n)) {
	    pred <- ilink(Xp %*% mrand[i, ])
	    opt[i] <- g$specs[which.max(pred)]
	  }
	  
      # mu <- mean(opt) # mean
      # cat("N*: ", mu, "\n")
      
      # se <- sd(opt) / sqrt(n) # standard error
      # cat("SE: ", se , "\n")
      
       quantile(opt, c(0.025, 0.975))
      
      # ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
      # cat("95% CI: ", ci , "\n")

		cat("\n P-spline smooth \n\n")
		gg  <- predict(HAC.ps, newdata = g, type = "response")
		Xp <- predict(HAC.ps, g, type = "lpmatrix") 
		beta <- coef(HAC.ps) 
		Vb <- vcov(HAC.ps)
		mrand <- mvrnorm(n, beta, Vb)
		opt <- rep(NA, n)
		ilink <- family(HAC.ps)$linkinv
		for (i in seq_len(n)) {
    	pred <- ilink(Xp %*% mrand[i, ])
    	opt[i] <- g$specs[which.max(pred)]
		}
		
        # mu <- mean(opt) # mean
        # cat("N*: ", mu, "\n")
        
        # se <- sd(opt) / sqrt(n) # standard error
        # cat("SE: ", se , "\n")
        
        quantile(opt, c(0.025, 0.975))
        
        # ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
        # cat("95% CI: ", ci , "\n")
		
		cat("\n Adaptive smooth \n\n")
		gg  <- predict(HAC.ad, newdata = g, type = "response")
		Xp <- predict(HAC.ad, g, type = "lpmatrix") 
		beta <- coef(HAC.ad) 
		Vb <- vcov(HAC.ad)
		mrand <- mvrnorm(n, beta, Vb)
		opt <- rep(NA, n)
		ilink <- family(HAC.ad)$linkinv
		for (i in seq_len(n)) {
		  pred <- ilink(Xp %*% mrand[i, ])
		  opt[i] <- g$specs[which.max(pred)]
		}

		
        # mu <- mean(opt) # mean
        # cat("N*: ", mu, "\n")
        
        # se <- sd(opt) / sqrt(n) # standard error
        # cat("SE: ", se , "\n")
        
        quantile(opt, c(0.025, 0.975))
        
        # ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
        # cat("95% CI: ", ci , "\n")

		}

	if (model == "SCAM") {
	  
	  cat("\n Monotonically increasing smooth \n\n")
	  gg  <- predict(HAC.mpi, newdata = g, type = "response") 
	  Xp <- predict(HAC.mpi, g, type = "lpmatrix") 
	  beta <- as.matrix(coef.scam(HAC.mpi)) 
	  Vb <- vcov.scam(HAC.mpi)
	  mrand <- mvrnorm(n, beta, Vb)
	  opt <- rep(NA, n)
	  ilink <- family(HAC.mpi)$linkinv
	  for (i in seq_len(n)) {
	    pred <- ilink(Xp %*% mrand[i, ])
	    opt[i] <- g$specs[which.max(pred)]
	  }
	  
	  mu <- mean(opt) 
	  cat("N*: ", mu, "\n")
      
      se <- sd(opt) / sqrt(n) # standard error
      cat("SE: ", se , "\n")
      
      ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
      cat("95% CI: ", ci , "\n")

	  cat("\n Concave smooth \n\n")
	  gg  <- predict(HAC.cv, newdata = g, type = "response") 
	  Xp <- predict(HAC.cv, g, type = "lpmatrix") 
	  beta <- as.matrix(coef.scam(HAC.cv)) 
	  Vb <- vcov.scam(HAC.cv)
	  mrand <- mvrnorm(n, beta, Vb)
	  opt <- rep(NA, n)
	  ilink <- family(HAC.cv)$linkinv
	  for (i in seq_len(n)) {
	    pred <- ilink(Xp %*% mrand[i, ])
	    opt[i] <- g$specs[which.max(pred)]
	  }

	  mu <- mean(opt) 
	  cat("N*: ", mu, "\n")
      
      se <- sd(opt) / sqrt(n) # standard error
      cat("SE: ", se , "\n")
      
      ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
      cat("95% CI: ", ci , "\n")
	
		cat("\n Monotonically increasing and concave smooth \n\n")
		gg  <- predict(HAC.micv, newdata = g, type = "response") 
		Xp <- predict(HAC.micv, g, type = "lpmatrix") 
		beta <- as.matrix(coef.scam(HAC.micv)) 
		Vb <- vcov.scam(HAC.micv)
		mrand <- mvrnorm(n, beta, Vb)
		opt <- rep(NA, n)
		ilink <- family(HAC.micv)$linkinv
		for (i in seq_len(n)) {
    			pred <- ilink(Xp %*% mrand[i, ])
    			opt[i] <- g$specs[which.max(pred)]
		}
		
		mu <- mean(opt) 
		cat("N*: ", mu, "\n")
        
        se <- sd(opt) / sqrt(n) # standard error
        cat("SE: ", se , "\n")
        
        ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
        cat("95% CI: ", ci , "\n")
        
		}

	if (model == "Krig") {
	
		cat("\n Matern covariance function \n\n")
		gg  <- predict(HAC.matern, newdata = g, type = "response") 
		Xp <- predict(HAC.matern, g, type = "lpmatrix") 
		beta <- coef(HAC.matern) 
		Vb <- vcov(HAC.matern)
		mrand <- mvrnorm(n, beta, Vb)
		opt <- rep(NA, n)
		ilink <- family(HAC.matern)$linkinv
		for (i in seq_len(n)) {
    			pred <- ilink(Xp %*% mrand[i, ])
    			opt[i] <- g$specs[which.max(pred)]
		}

	
		mu <- mean(opt) 
		cat("N*: ", mu, "\n")
        
        se <- sd(opt) / sqrt(n) # standard error
        cat("SE: ", se , "\n")
        
        ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
        cat("95% CI: ", ci , "\n")
		
		cat("\n Spherical covariance function \n\n")
		gg  <- predict(HAC.sph, newdata = g, type = "response") 
		Xp <- predict(HAC.sph, g, type = "lpmatrix") 
		beta <- coef(HAC.sph) 
		Vb <- vcov(HAC.sph)
		mrand <- mvrnorm(n, beta, Vb)
		opt <- rep(NA, n)
		ilink <- family(HAC.sph)$linkinv
		for (i in seq_len(n)) {
    			pred <- ilink(Xp %*% mrand[i, ])
    			opt[i] <- g$specs[which.max(pred)]
		}
	
		mu <- mean(opt) 
		cat("N*: ", mu, "\n")
        
        se <- sd(opt) / sqrt(n) # standard error
        cat("SE: ", se , "\n")
        
        ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
        cat("95% CI: ", ci , "\n")
		
		cat("\n Exponential covariance function \n\n")
		gg  <- predict(HAC.exp, newdata = g, type = "response") 
		Xp <- predict(HAC.exp, g, type = "lpmatrix") 
		beta <- coef(HAC.exp) 
		Vb <- vcov(HAC.exp) 
		mrand <- mvrnorm(n, beta, Vb) 
		opt <- rep(NA, n)
		ilink <- family(HAC.exp)$linkinv
		for (i in seq_len(n)) {
		  pred <- ilink(Xp %*% mrand[i, ])
		  opt[i] <- g$specs[which.max(pred)]
		}

		mu <- mean(opt) 
		cat("N*: ", mu, "\n")
        
        se <- sd(opt) / sqrt(n) # standard error
        cat("SE: ", se , "\n")
        
        ci <- mu + c(-1, 1) * qnorm(0.975) * se # confidence interval
        cat("95% CI: ", ci , "\n")
 
		
	}

}
