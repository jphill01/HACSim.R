HAC.simbestaic <- function(model) {
  
  # Find the best single model based on the category searching for
  
  if (model == "GAM") {
    minAIC <- c("tp", "cr", "ps","ad")[which.min(c(HAC.tp$aic, HAC.cr$aic, HAC.ps$aic, HAC.ad$aic))]
  } else if (model == "SCAM") {
    minAIC <- c("mpi", "cv", "micv")[which.min(c(HAC.mpi$aic, HAC.cv$aic, HAC.micv$aic))]
  } else if (model == "Krig") {
    minAIC <- c("matern", "sph", "exp")[which.min(c(HAC.matern$aic, HAC.sph$aic, HAC.exp$aic))]
  } else if (model == "All") {
    minAIC <- c("tp", "cr", "ps","ad",
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
