\name{HAC.simrep}
\alias{HAC.simrep}

\title{Run a simulation}

\description{Runs the HACSim algorithm to iteratively extrapolate species' haplotype 
accumulation curves 
}

\usage{HAC.simrep(filename)}

\arguments{\item{filename}{Name of file where simulation results are saved}
}

\examples{
  ## Simulate hypothetical species ##
  
  N <- 100 # total number of sampled individuals
  Hstar <- 10 # total number of haplotypes
  probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution
  
  HACSObj <- HACHypothetical(N, Hstar, probs, filename = "output") # outputs a CSV file called "output.csv"
  
  ## Simulate hypothetical species - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 1000, p = 0.95, subsample = TRUE, prop = 0.25, conf.level = 0.95, filename = "output")
  
  ## Simulate hypothetical species and all paramaters changed - subsampling ##
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.95, filename = "output")
  
  HAC.simrep(HACSObj) # runs a simulation
  
  
  ## Simulate real species ##
  
  \dontrun{
    ## Simulate real species ##
    HACSObj <- HACReal(filename = "output") # outputs file called "output.csv"
    
    ## Simulate real species - subsampling ##
    HACSObj <- HACReal(subsample = TRUE, prop = 0.15, conf.level = 0.95, filename = "output")
    
    ## Simulate real species and all parameters changed - subsampling ##
    HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, conf.level = 0.99, filename = "output")
    
    HAC.simrep(HACSObj) # user prompted to select appropriate FASTA file of aligned/trimmed DNA sequences prior to running a simulation}