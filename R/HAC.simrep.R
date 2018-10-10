### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function(filename = "output") {
  ptm <- proc.time()
  assign("iters", 1, .GlobalEnv)
  
  df <- data.frame(matrix(ncol = 10, nrow = 0))
  x <- c("Mean number of haplotypes sampled",
         "Lower 95% confidence limit for number of haplotypes recovered",
         "Upper 95% confidence limit for number of haplotypes recovered",
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
               p = p,
               subset.haps = subset.haps,
               prop.haps = prop.haps,
               subset.seqs = subset.seqs,
               prop.seqs = prop.seqs,
               input.seqs = input.seqs,
               sim.seqs = sim.seqs,
               num.seqs = num.seqs,
               length.seqs = length.seqs,
               subst.model = subst.model,
               nucl.freq = nucl.freq,
               mu.rate = mu.rate,
               transi.rate = transi.rate,
               transv.rate = transv.rate,
               num.pts = 10,
               prop.pts = NULL,
               df = df
  )
  
  amt <- proc.time() - ptm
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("\n \n \n Desired level of H* has not yet been reached \n")
  } else {
    cat("\n \n \n Desired level of H* has been reached. \n \n The algorithm converged after", 
        iters, "iterations and took", amt[3], "s.", "\n")
  }
  
  while (R < p) {
    df <- HAC.sim(N = ceiling(Nstar), 
                  Hstar = Hstar, 
                  probs = probs,
                  perms = perms,
                  p = p,
                  subset.haps = subset.haps,
                  prop.haps = prop.haps,
                  subset.seqs = subset.seqs,
                  prop.seqs = prop.seqs,
                  num.pts = 10,
                  prop.pts = NULL,
                  df = df)
    
    assign("iters", iters + 1, .GlobalEnv)
    amt <- proc.time() - ptm
    
    ## Check whether desired level of haplotype recovery has been reached ##
    
    if (R < p) {
      cat("\n \n \n Desired level of H* has not yet been reached \n")
    } else {
      cat("\n \n Desired level of H* has been reached. \n \n ---------- Finished. ----------
          \n The algorithm converged after", iters, "iterations and took", amt[3], "s.", 
          "\n \n The estimate of sampling sufficiency for", p * 100,"% haplotype recovery is:" , max(d$specs))
    }
    
  }
  
  write.csv(df, file = paste(filename, ".csv", sep = ""))
} # end HAC.simrep
