### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function(filename = "output") {
  ptm <- proc.time()
  assign("iters", 1, .GlobalEnv)
  
  df <- data.frame(matrix(ncol = 8, nrow = 0))
  x <- c("Mean number of haplotypes sampled", 
         "Mean number of haplotypes not sampled", 
         "Proportion of haplotypes (specimens) sampled", 
         "Proportion of haplotypes (specimens) not sampled",
         "Mean value of N*", 
         "Mean number of specimens not sampled",
         "Haplotype accumulation curve slope", 
         "Mean number of specimens required to observe one new haplotype")
  colnames(df) <- x
  
  df <- HAC.sim(N = N, 
               Hstar = Hstar, 
               probs = probs, 
               perms = perms,
               K = K,
               p = p,
               subset.haps = subset.haps,
               prop.haps = prop.haps,
               input.seqs = input.seqs,
               subset.seqs = subset.seqs,
               prop.seqs = prop.seqs,
               num.pts = 10,
               prop.pts = NULL,
               df = df
  )
  
  amt <- proc.time() - ptm
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("Desired level of H* has not yet been reached \n")
  } else {
    cat("\n Desired level of H* has been reached. \n \n The algorithm converged after", 
        iters, "iterations and took", amt[3], "s.", "\n")
  }
  
  while (R < p) {
    df <- HAC.sim(N = ceiling(Nstar), 
                  Hstar = Hstar, 
                  probs = probs,
                  perms = perms,
                  K = K,
                  p = p,
                  subset.haps = subset.haps,
                  prop.haps = prop.haps,
                  subset.seqs = subset.seqs,
                  prop.seqs = prop.seqs,
                  num.pts = 10,
                  prop.pts = NULL,
                  df = df
    )
    
    assign("iters", iters + 1, .GlobalEnv)
    amt <- proc.time() - ptm
    
    ## Check whether desired level of haplotype recovery has been reached ##
    
    if (R < p) {
      cat("Desired level of H* has not yet been reached \n")
    } else {
      cat("\n Desired level of H* has been reached. \n \n The algorithm converged after", 
          iters + 1, "iterations and took", amt[3], "s.", "\n")
    }
    
    
  }
  
  write.csv(df, file = paste(filename, ".csv", sep = ""))
}
