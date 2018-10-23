### Semiparametric Model Fitting ###

HAC.simaic <- function(model = c("GAM", "SCAM", "Krig", "Best", "All"), k = 10) {
  
  if (model == "All") {
    
    cat("\n Thin plate smooth (tp) \n")
    cat("\n AIC: " ,  HAC.tp$aic, "\n")
    
    cat("\n Cubic spline smooth (cr)\n")
    cat("\n AIC: " ,  HAC.cr$aic, "\n")
    
    cat("\n P-spline smooth (ps) \n")
    cat("\n AIC: " ,  HAC.ps$aic, "\n")
    
    cat("\n Adaptive smooth (ad) \n")
    cat("\n AIC: " ,  HAC.ad$aic, "\n")
    
    cat("\n Monotonically increasing smooth (mpi) \n")	
    cat("\n AIC: " ,  HAC.mpi$aic, "\n")
    
    cat("\n Concave smooth (cv)\n")	
    cat("\n AIC: " ,  HAC.cv$aic, "\n")
    
    cat("\n Monotonically increasing and concave smooth (micv) \n")	
    cat("\n AIC: " ,  HAC.micv$aic, "\n")
    
    cat("\n Matern covariance function (Matern) \n")
    cat("\n AIC: " ,  HAC.matern$aic, "\n")
    
    cat("\n Spherical covariance function (sph) \n")
    cat("\n AIC: " ,  HAC.sph$aic, "\n")
    
    cat("\n Exponential covariance function (exp) \n")
    cat("\n AIC: " ,  HAC.exp$aic, "\n")
    
  }
  
  if ((model == "GAM") || (model == "SCAM") || (model == "Krig") || (model == "Best")) {
  
  # Find the best single model based on the category searching for
  
  if (model == "GAM") {
    minAIC <- c("tp", "cr", "ps","ad")[which.min(c(HAC.tp$aic, HAC.cr$aic, HAC.ps$aic, HAC.ad$aic))]
  } else if (model == "SCAM") {
    minAIC <- c("mpi", "cv", "micv")[which.min(c(HAC.mpi$aic, HAC.cv$aic, HAC.micv$aic))]
  } else if (model == "Krig") {
    minAIC <- c("matern", "sph", "exp")[which.min(c(HAC.matern$aic, HAC.sph$aic, HAC.exp$aic))]
  } else if (model == "Best") {
    minAIC <- c("tp", "cr", "ps", "ad",
                "mpi", "cv", "micv",
                "matern","sph", "exp")[which.min(c(HAC.tp$aic, HAC.cr$aic, HAC.ps$aic, HAC.ad$aic,
                                                   HAC.mpi$aic, HAC.cv$aic, HAC.micv$aic,
                                                   HAC.matern$aic, HAC.sph$aic, HAC.exp$aic))]
  }
  
  # Print the name of the model chosen, and its corresponding AIC
  
  switch(minAIC,
         tp = {cat("\n Thin plate spline smooth \n\n AIC: ", HAC.tp$aic)},
         cr = {cat("\n Cubic smooth \n\n AIC: ", HAC.cr$aic)},
         ps = {cat("\n P-spline smooth \n\n AIC: ", HAC.ps$aic)},
         ad = {cat("\n Adaptive smooth \n\n AIC: ", HAC.ad$aic)},
         mpi = {cat("\n Monotonically increasing smooth \n\n AIC: ", HAC.mpi$aic)},
         cv = {cat("\n Concave smooth \n\n AIC: ", HAC.cv$aic)},
         micv = {cat("\n Monotonically increasing and concave smooth \n\n AIC: ", HAC.micv$aic)},
         matern = {cat("\n Matern covariance function \n\n AIC: ", HAC.matern$aic)},
         sph = {cat("\n Spherical covariance function \n\n AIC: ", HAC.sph$aic)},
         exp = {cat("\n Exponential covariance function \n\n AIC: ", HAC.exp$aic)}
  )
    
  }

}
