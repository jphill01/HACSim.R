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
  
  ## Calculate slope of final curve using last n points (or proportion of points) on curve
  ## perms must be large enough to ensure monotonicity and a non-negative slope
  
  if ((!is.null(prop.pts)) && (is.null(num.pts))) { 
    lin.reg <- lm(means ~ specs, data = tail(d, n = ceiling(prop.pts * nrow(d))))
  }
  
  if ((!is.null(num.pts)) && (is.null(prop.pts))) {
    lin.reg <- lm(means ~ specs, data = tail(d, n = num.pts))
  }
  
  beta1 <- abs(coef(lin.reg)[[2]])
  
  out <- calibrate(lin.reg, y0 = R * Hstar, interval = "Wald", mean.response = FALSE)
  
  ## Check whether desired level of haplotype recovery has been reached ##
  
  if (R < p) {
    cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
  } else {
    cat("\n \n \n Desired level of haplotype recovery has been reached \n \n \n ---------- Finished. ----------
        \n --- Summary Statistics Results ---
        \n The initial guess for sampling sufficiency was N =", paste0(N, "."),
        "\n \n The algorithm converged after", iters, "iterations and took", amt[3], "s.", 
        "\n \n The estimate of sampling sufficiency for p =", paste0(p * 100, "%"), "haplotype recovery is N* = ", max(d$specs), "individuals.",
        "\n \n The number of additional specimens required to be sampled for p = ", paste0(p * 100, "%"), "haplotype recovery is \n N* - N =", max(d$specs) - N, "individuals.",
        "\n \n --- Linear Model Results --- \n \n", 
        "Haplotype accumulation curve slope: ", beta1,
        "\n Mean number of specimens required to observe one new haplotype: ", 1 / beta1,
        "\n \n Mean value of N*: ", out$estimate, "( 95% CI:", paste(ceiling(out$lower), ceiling(out$upper), sep = "-"), ")")
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
    
    ## Calculate slope of final curve using last n points (or proportion of points) on curve
    ## perms must be large enough to ensure monotonicity and a non-negative slope
    
    if ((!is.null(prop.pts)) && (is.null(num.pts))) { 
      lin.reg <- lm(means ~ specs, data = tail(d, n = ceiling(prop.pts * nrow(d))))
    }
    
    if ((!is.null(num.pts)) && (is.null(prop.pts))) {
      lin.reg <- lm(means ~ specs, data = tail(d, n = num.pts))
    }
    
    beta1 <- abs(coef(lin.reg)[[2]])
    
    out <- calibrate(lin.reg, y0 = R * Hstar, interval = "Wald", mean.response = FALSE)
    
    ## Check whether desired level of haplotype recovery has been reached ##
    
    if (R < p) {
      cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
    } else {
      cat("\n \n \n Desired level of haplotype recovery has been reached \n \n \n ---------- Finished. ----------
          \n --- Summary Statistics Results ---
          \n The initial guess for sampling sufficiency was N = ", paste0(N, "."),
          "\n \n The algorithm converged after", iters, "iterations and took", amt[3], "s.", 
          "\n \n The estimate of sampling sufficiency for p =", paste0(p * 100, "%"), "haplotype recovery is N* = ", max(d$specs), "individuals.",
          "\n \n The number of additional specimens required to be sampled for p =", paste0(p * 100, "%"), "haplotype recovery is \n N* - N = ",  max(d$specs) - N, "individuals.",
          "\n \n --- Linear Model Results --- \n \n", 
          "Haplotype accumulation curve slope: ", beta1,
          "\n Mean number of specimens required to observe one new haplotype: ", 1 / beta1,
          "\n \n Mean value of N*: ", out$estimate, "( 95% CI:", paste(ceiling(out$lower), ceiling(out$upper), sep = "-"), ")")
    }
    
  }
  
  write.csv(df, file = paste(filename, ".csv", sep = ""))
} # end HAC.simrep
