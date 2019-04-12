envr <- NULL 

.onLoad <- function(...) {
  envr <<- new.env()  # when package is loaded, create new environment to store needed variables 
}