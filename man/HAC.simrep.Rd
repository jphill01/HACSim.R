\name{HAC.simrep}
\alias{HAC.simrep}

\title{Run a simulation of haplotype accumulation curves}

\description{Runs the HACSim algorithm to iteratively extrapolate haplotype 
accumulation curves for hypothetical or real species
}

\usage{HAC.simrep(HACSObj)}

\arguments{\item{HACSObj}{Name of variable where simulation results are stored}
}

\note{When simulating real species via HACReal(...), a pop-up window will appear prompting the user to select an intraspecific FASTA file of aligned/trimmed DNA sequences. The alignment must not contain missing or ambiguous nucleotides; i.e., it should only contain A, C, G or T.

HACSObj can be altered to any name desired by the user.

}

\examples{
  ## Simulate hypothetical species ##
  
  N <- 100 # total number of sampled individuals
  Hstar <- 10 # total number of haplotypes
  probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution
  
  HACSObj <- HACHypothetical(N = N, Hstar = Hstar , probs = probs, filename = "output") # outputs a CSV file called "output.csv"
  
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
    
    HAC.simrep(HACSObj) # user prompted to select appropriate FASTA file of aligned/trimmed DNA sequences prior to running a simulation
    }
}