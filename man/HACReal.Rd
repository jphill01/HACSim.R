\name{HACReal}
\alias{HACReal}

\title{Set up object to simulate real species}

\usage{HACReal(...)}

\examples{
  \dontrun{
    ## Simulate real species ##
    HACSObj <- HACReal(filename = "output") # outputs file called "output.csv"
    
    ## Simulate real species - subsampling ##
    HACSObj <- HACReal(subsample = TRUE, prop = 0.15, conf.level = 0.95, filename = "output")
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.99, filename = "output")
}
