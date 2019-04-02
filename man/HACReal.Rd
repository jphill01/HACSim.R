\name{HACReal}
\alias{HACReal}

\title{Set up an object to simulate haplotype accumulation for a real species}

\description{Helper function which creates an object containing necessary information to run a simulation of haplotype accumulation for a real species of interest}

\usage{HACReal(...)}

\arguments{\item{p}{Proportion of haplotypes to recover}
\item{perms}{Number of permutations (replications)}
\item{subsample}{Is a subsample of DNA sequences desired (TRUE/FALSE)?}
\item{prop}{If subsample = TRUE, the proportion of DNA sequences to subsample}
\item{conf.level}{Desired confidence level for gaphical output and interval estimation}
\item{filename}{Name of file where simulation results are to be saved}
}

\examples{
  \dontrun{
    ## Simulate real species ##
    HACSObj <- HACReal(filename = "output") # outputs file called "output.csv"
    
    ## Simulate real species - subsampling ##
    HACSObj <- HACReal(subsample = TRUE, prop = 0.15, conf.level = 0.95, filename = "output")
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.99, filename = "output")
}
