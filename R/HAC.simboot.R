### Bootstrap Simulation ###

HAC.simboot <- function(model = c("GAM", "SCAM", "Krig"), k = 10, bootType="Bisect") {
    "Finds the name for the best model"
    HAC.findMinAIC <- function() {
        maxAIC <- c("tp", "cr", "ps","ad",
                    "mpi", "cv", "micv",
                    "matern","sph", "exp")[which.min(c(HAC.tp$aic,HAC.cr$aic,HAC.ps$aic,HAC.ad$aic,
                                                       HAC.mpi$aic,HAC.cv$aic,HAC.micv$aic,
                                                       HAC.matern$aic,HAC.sph$aic,HAC.exp$aic))]
    }
    # Sets the models to be tested, alternatively can set to a specific model with the shortened model name as the model argument
    if (model == "GAM") {
      model <- c('tp','cr','ps','ad')
    } else if (model == "SCAM") { 
      model <- c('mpi','cv','micv')
    } else if (model == "Krig") {
      model <- c('matern','sph','exp')
    } else if (model == "best") { 
      model <- HAC.findMinAIC()
    } else if (model == "all") { 
      model <- c("tp", "cr", "ps", "ad", "mpi", "cv", "micv", "matern","sph", "exp")
    }

    for (i in model) {
        # Set up the paramaters for the bootstrapping based on the model currently being tested
        if (i == 'tp') {
          modelName = "tp"; currModel = HAC.tp; displayName = "Thin plate smooth (tp)"; modelType = "GAM"
        } else if(i == 'cr') {
          modelName = "cr"; currModel = HAC.cr; displayName = "Cubic spline smooth (cr)"; modelType = "GAM"
        } else if (i == 'ps') {
          modelName = "ps"; currModel = HAC.ps; displayName = "P-spline smooth (ps)"; modelType = "GAM"
        } else if (i == 'ad') {
          modelName = "ad"; currModel = HAC.ad; displayName = "Adaptive smooth (ad)"; modelType = "GAM"
        } else if (i == 'mpi') {
          modelName = "mpi"; currModel = HAC.mpi; displayName = "Monotonically increasing smooth (mpi)"; modelType = "SCAM"
        } else if (i == 'cv') {
          modelName = "cv"; currModel = HAC.cv; displayName = "Concave smooth (cv)"; modelType = "SCAM"
        } else if (i == 'micv') {
          modelName = "micv"; currModel = HAC.micv; displayName = "Monotonically increasing and concave smooth (micv)"; modelType = "SCAM"
        ###The modelType is GAM because it has similar arguments to the other GAM models, this is just for convenience (no mVal for it)
        } else if (i == 'matern') {
          modelName = "gp"; currModel = HAC.matern; displayName = "Matern covariance function"; modelType = "GAM"
        } else if (i == 'sph') {
          modelName = "gp"; currModel = HAC.sph; displayName = "Spherical covariance function (sph)"; modelType = "Krig"; mVal = 1
        } else if (i == 'exp') {
          modelName = "gp"; currModel = HAC.exp; displayName = "Exponential covariance function (exp)"; modelType = "Krig"; mVal = 2
        } else {
          stop("Invalid model")
          }
        
        cat("\n" ,displayName, "\n")
        res <- resid(currModel) - mean(resid(currModel)) # centre the residuals
        n <- length(res)
        boot.data <- data.frame(d, res = res, fit = fitted(currModel))
        boot.fun <- function(data, i) {
            ## Differenced based on the group of models. 
            if (modelType == "GAM"){
              boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = modelName, k = k), method = "GCV.Cp", optimizer = c("outer", "bfgs"), data = data)
            } else if (modelType == "SCAM"){
              boot.fit <- scam(boot.data$means + res[i] ~ s(specs, bs = modelName, k = k), data = data)
            } else if (modelType =="Krig") {
                boot.fit <- gam(boot.data$means + res[i] ~ s(specs, bs = modelName, k = k, m = mVal), method = "GCV.Cp", optimizer = c("outer", "bfgs"), data = data)
            }
            
            # Simulate the correct variance
            Y0 <- R * Hstar + sample(data$res, size = 1, replace = TRUE)
            # Make sure the original estimate also gets returned
            
            #The only lines that differ between the newton and bisection methods
            if (bootType == "Bisect") {
                if (all(i == 1:n)) {
                    inv.predict(currModel, y = R*Hstar, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
                } else {
                    inv.predict(boot.fit, y = Y0, x.name = "specs", lower = 1, upper = ceiling(Nstar), interval = FALSE)[1L]
                }
            }
            else if (bootType == "Newton") {
                if (all(i == 1:n)) {
                    inv.predict.newton(currModel, y = R*Hstar, x.name = "specs", start = N, interval = FALSE)
                } else {
                    inv.predict.newton(boot.fit, y = Y0, x.name = "specs", start = N, interval = FALSE)
                }
            }
            else {
                stop("Not a valid bootstrap model type")
            }
        }
        
        res <- boot(boot.data, boot.fun, R = 1000)
        print(res)
        rel.bias <- (mean(res$t, na.rm = TRUE) - res$t0) / res$t0 # relative bias
        cat("\n Percent relative bias: ", rel.bias*100) 
        plot(res)
        cat("\n\n")
        print(boot.ci(res, type = "all"))
  }
}
