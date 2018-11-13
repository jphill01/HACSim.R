### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function(filename = "output") {
  ptm <- proc.time()
  assign("iters", 1, .GlobalEnv)
  
  df <- data.frame(matrix(ncol = 8, nrow = 0))
  x <- c("Mean number of haplotypes sampled",
         "Lower 95% confidence limit for number of haplotypes recovered",
         "Upper 95% confidence limit for number of haplotypes recovered",
         "Mean number of haplotypes not sampled", 
         "Proportion of haplotypes (specimens) sampled", 
         "Proportion of haplotypes (specimens) not sampled",
         "Mean value of N*", 
         "Mean number of specimens not sampled")
  colnames(df) <- x
  
  cat("\n Simulating haplotype accumulation...")
  
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
               subst.model = subst.model,
               df = df
  )
  
  amt <- proc.time() - ptm
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
  } else {
    cat("\n \n \n Desired level of haplotype recovery has been reached \n \n \n ---------- Finished. ----------
        \n The initial guess for sampling sufficiency was N = ", paste0(N, "."),
        "\n \n The algorithm converged after", iters, "iterations and took", amt[3], "s.", 
        "\n \n The estimate of sampling sufficiency for p =", paste0(p * 100, "%"), "haplotype recovery is N* = ", max(d$specs), "individuals.",
        "\n \n The number of additional specimens required to be sampled for p =", paste0(p * 100, "%"), "haplotype recovery is \n N* - N = ",  max(d$specs) - N, "individuals.")
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
                  df = df)
    
    assign("iters", iters + 1, .GlobalEnv)
    amt <- proc.time() - ptm
    
    ## Check whether desired level of haplotype recovery has been reached ##
    
    if (R < p) {
      cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
    } else {
      cat("\n \n \n Desired level of haplotype recovery has been reached \n \n \n ---------- Finished. ----------
          \n The initial guess for sampling sufficiency was N = ", paste0(N, "."),
          "\n \n The algorithm converged after", iters, "iterations and took", amt[3], "s.", 
          "\n \n The estimate of sampling sufficiency for p =", paste0(p * 100, "%"), "haplotype recovery is N* = ", max(d$specs), "individuals.",
          "\n \n The number of additional specimens required to be sampled for p =", paste0(p * 100, "%"), "haplotype recovery is \n N* - N = ",  max(d$specs) - N, "individuals.")
    }
    
  }
  
  write.csv(df, file = paste(filename, ".csv", sep = ""))
} # end HAC.simrep
