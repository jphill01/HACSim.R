### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function(HACSObject) {
    assign("N", HACSObject$N, envir = .GlobalEnv)
    assign("Hstar", HACSObject$Hstar, envir = .GlobalEnv)
    assign("probs", HACSObject$probs, envir = .GlobalEnv)
    perms <- HACSObject$perms
    p <- HACSObject$p
    conf.level <- HACSObject$conf.level
    subset.haps <- HACSObject$subset.haps
    prop.haps <- HACSObject$prop.haps
    subset.seqs <- HACSObject$subset.seqs
    prop.seqs <- HACSObject$prop.seqs
    input.seqs <- HACSObject$input.seqs
    filename <- HACSObject$filename
  
  assign("ptm", proc.time(), envir = .GlobalEnv)
  assign("iters", 1, .GlobalEnv)
  
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  x <- c("Mean number of haplotypes sampled",
         "Mean number of haplotypes not sampled", 
         "Proportion of haplotypes sampled", 
         "Proportion of haplotypes not sampled",
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
               conf.level = conf.level,
               df = df
  )
  
  amt <- proc.time() - ptm
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
  } else {
    cat("\n \n \n Desired level of haplotype recovery has been reached \n \n \n ---------- Finished. ----------
        \n The initial guess for sampling sufficiency was N = ", paste0(N), "individuals",
        "\n \n The algorithm converged after", iters, "iterations and took", amt[3], "s", 
        "\n \n The estimate of sampling sufficiency for p =", paste0(p * 100, "%"), "haplotype recovery is N* = ", (low + high) / 2, "individuals (", paste0(conf.level * 100, "%"), "CI:", paste(low, high, sep = "-"), ")",
        "\n \n The number of additional specimens required to be sampled for p =", paste0(p * 100, "%"), "haplotype recovery is \n N* - N = ",  ((low + high) / 2) - N, "individuals")
  }
  
  while (R < p) {
    df <- HAC.sim(N = Nstar, 
                  Hstar = Hstar, 
                  probs = probs,
                  perms = perms,
                  p = p,
                  subset.haps = subset.haps,
                  prop.haps = prop.haps,
                  subset.seqs = subset.seqs,
                  prop.seqs = prop.seqs,
                  conf.level = conf.level,
                  df = df)
    
    assign("iters", iters + 1, .GlobalEnv)
    amt <- proc.time() - ptm
    
    ## Check whether desired level of haplotype recovery has been reached ##
    
    if (R < p) {
      cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
    } else {
      cat("\n \n \n Desired level of haplotype recovery has been reached \n \n \n ---------- Finished. ----------
          \n The initial guess for sampling sufficiency was N = ", paste0(N), "individuals",
          "\n \n The algorithm converged after", iters, "iterations and took", amt[3], "s", 
          "\n \n The estimate of sampling sufficiency for p =", paste0(p * 100, "%"), "haplotype recovery is N* = ", (low + high) / 2, "individuals (",  paste0(conf.level * 100, "%"), "CI:", paste(low, high, sep = "-"), ")",
          "\n \n The number of additional specimens required to be sampled for p =", paste0(p * 100, "%"), "haplotype recovery is \n N* - N = ", ((low + high) / 2) - N, "individuals")
    }
    
  } 
  if (!is.null(filename))
    write.csv(df, file = paste(filename, ".csv", sep = ""))
} # end HAC.simrep
