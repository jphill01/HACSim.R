HAC.simmodels <- function(k = 10) {
  
  ## GAMs ##

  assign("HAC.tp", gam(means ~ s(specs, bs = "tp", k = k), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  assign("HAC.cr", gam(means ~ s(specs, bs = "cr", k = k), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  assign("HAC.ps", gam(means ~ s(specs, bs = "ps", k = k), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  assign("HAC.ad", gam(means ~ s(specs, bs = "ad", k = k), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  
  ## SCAMs ##
  
  assign("HAC.mpi", scam(means ~ s(specs, bs = "mpi", k = k), data = d), envir = .GlobalEnv)
  assign("HAC.cv", scam(means ~ s(specs, bs = "cv", k = k), data = d), envir = .GlobalEnv)
  assign("HAC.micv", scam(means ~ s(specs, bs = "micv", k = k), data = d), envir = .GlobalEnv)
  
  ## Kriging ##
  
  assign("HAC.matern", gam(means ~ s(specs, bs = "gp", k = k), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  assign("HAC.sph", gam(means ~ s(specs, bs = "gp", k = k, m = 1), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  assign("HAC.exp", gam(means ~ s(specs, bs = "gp", k = k, m = 2), optimizer = c("outer", "bfgs"), data = d), envir = .GlobalEnv)
  
}