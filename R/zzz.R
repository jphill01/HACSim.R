envr <- NULL 

.onLoad <- function(...) {
  envr <<- new.env()  # when package is loaded, create new environment to store needed variables 
}

.onAttach <- function(...) {
  packageStartupMessage("This is HACSim 1.0.0")
}