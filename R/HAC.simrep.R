### HAC Simulation Iteration ###

## Run HAC Simulator until convergence (saturation) is reached ##

HAC.simrep <- function(HACSObject) {
  assign("N", HACSObject$N, envir = envr)
  assign("Hstar", HACSObject$Hstar, envir = envr)
  assign("probs", HACSObject$probs, envir = envr)
  assign("perms", HACSObject$perms, envir = envr)
  assign("p", HACSObject$p, envir = envr)
  assign("conf.level", HACSObject$conf.level, envir = envr)
  assign("ci.type", HACSObject$ci.type, envir = envr)
  assign("subset.haps", HACSObject$subset.haps, envir = envr)
  assign("prop.haps", HACSObject$prop.haps, envir = envr)
  assign("subset.seqs", HACSObject$subset.seqs, envir = envr)
  assign("prop.seqs", HACSObject$prop.seqs, envir = envr)
  assign("input.seqs", HACSObject$input.seqs, envir = envr)
  assign("progress", HACSObject$progress, envir = envr)
  assign("num.iters", HACSObject$num.iters, envir = envr)
  assign("filename", HACSObject$filename, envir = envr)

  assign("ptm", proc.time(), envir = envr)
  assign("iters", 1, envir = envr)

  df <- data.frame(matrix(ncol = 6, nrow = 0))
  x <- c(
    "Mean number of haplotypes sampled",
    "Mean number of haplotypes not sampled",
    "Proportion of haplotypes sampled",
    "Proportion of haplotypes not sampled",
    "Mean value of N*",
    "Mean number of specimens not sampled"
  )
  colnames(df) <- x

  cat("\n Simulating haplotype accumulation...")

  df <- HAC.sim(
    N = envr$N,
    Hstar = envr$Hstar,
    probs = envr$probs,
    perms = envr$perms,
    p = envr$p,
    subset.haps = envr$subset.haps,
    prop.haps = envr$prop.haps,
    subset.seqs = envr$subset.seqs,
    prop.seqs = envr$prop.seqs,
    input.seqs = envr$input.seqs,
    conf.level = envr$conf.level,
    ci.type = envr$ci.type,
    progress = envr$progress,
    num.iters = envr$num.iters,
    df = df
  )

  amt <- proc.time() - envr$ptm

  ## Check whether desired level of haplotype recovery has been reached ##

  if (envr$progress == TRUE) {
    if (envr$R < envr$p) {
      cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
    } else {
      cat(
        "\n \n \n Desired level of haplotype recovery has been reached 
        \n \n \n ---------- Finished. ----------
        \n The initial guess for sampling sufficiency was N = ", envr$N, "individuals represented by \n H* = ", envr$Hstar, "haplotypes",
        "\n \n The fraction of species' haplotype diversity captured from sampling N = ", envr$N, "individuals is R = ", envr$df.out[1, 3],
        "\n \n The algorithm converged after", envr$iters, "iterations with", envr$perms, "permutations and took", amt[3], "s",
        "\n \n The estimate of sampling sufficiency for p =", paste0(envr$p * 100, "%"), "haplotype recovery is N* = ", envr$Nstar - envr$X, "individuals \n (", paste0(envr$conf.level * 100, "%"), "CI:", paste(envr$Nstar.low, envr$Nstar.high, sep = "-"), ")",
        "\n \n The estimated fraction of species' haplotype diversity captured from sampling N* = ", envr$Nstar - envr$X, "individuals is R =", signif(envr$R, digits = 3), "(", paste0(envr$conf.level * 100, "%"), "CI:", paste(signif(tail(envr$R.low, n = 1), digits = 3), signif(tail(envr$R.high, n = 1), digits = 3), sep = "-"), ")",
        "\n \n The number of additional specimens required to be sampled for p =", paste0(envr$p * 100, "%"), "haplotype recovery is N* - N = ", envr$Nstar - envr$X - envr$N, "individuals (", paste0(envr$conf.level * 100, "%"), "CI:", paste(envr$Nstar.low - envr$N, envr$Nstar.high - envr$N, sep = "-"), ") 
        \n \n -------------------------------"
      )
    }
  }
  
  if (is.null(envr$num.iters)) {
    while (envr$R < envr$p) {
      df <- HAC.sim(
        N = ceiling(envr$Nstar), # round sample size up
        Hstar = envr$Hstar,
        probs = envr$probs,
        perms = envr$perms,
        p = envr$p,
        subset.haps = envr$subset.haps,
        prop.haps = envr$prop.haps,
        subset.seqs = envr$subset.seqs,
        prop.seqs = envr$prop.seqs,
        conf.level = envr$conf.level,
        ci.type = envr$ci.type,
        num.iters = envr$num.iters,
        progress = envr$progress,
        df = df
      )

      assign("iters", envr$iters + 1, envr)
      amt <- proc.time() - envr$ptm

      ## Check whether desired level of haplotype recovery has been reached ##

      if (envr$progress == TRUE) {
        if (envr$R < envr$p) {
          cat("\n \n \n Desired level of haplotype recovery has not yet been reached \n")
        } else {
          cat(
            "\n \n \n Desired level of haplotype recovery has been reached 
            \n \n \n ---------- Finished. ----------
            \n The initial guess for sampling sufficiency was N = ", envr$N, "individuals represented by H* = ", envr$Hstar, "haplotypes",
            "\n \n The fraction of species' haplotype diversity captured from sampling N = ", envr$N, "individuals is R = ", envr$df.out[1, 3],
            "\n \n The algorithm converged after", envr$iters, "iterations with", envr$perms, "permutations and took", amt[3], "s",
            "\n \n The estimate of sampling sufficiency for p =", paste0(envr$p * 100, "%"), "haplotype recovery is N* = ", envr$Nstar - envr$X, "individuals \n (", paste0(envr$conf.level * 100, "%"), "CI:", paste(envr$Nstar.low, envr$Nstar.high, sep = "-"), ")",
            "\n \n The estimated fraction of species' haplotype diversity captured from sampling N* = ", envr$Nstar - envr$X, "individuals is R =", signif(envr$R, digits = 3), "(", paste0(envr$conf.level * 100, "%"), "CI:", paste(signif(tail(envr$R.low, n = 1), digits = 3), signif(tail(envr$R.high, n = 1), digits = 3), sep = "-"), ")",
            "\n \n The number of additional specimens required to be sampled for p =", paste0(envr$p * 100, "%"), "haplotype recovery is N* - N = ", envr$Nstar - envr$X - envr$N, "individuals (", paste0(envr$conf.level * 100, "%"), "CI:", paste(envr$Nstar.low - envr$N, envr$Nstar.high - envr$N, sep = "-"), ") 
            \n \n -------------------------------"
          )
        }
      }
    }
  }
  
  assign("df.out", df, envir = envr)
  
  cat("\n \n")

  message("\n \n Type envr$ to extract simulation parameters of interest (see documentation for details)")

    if (!is.null(envr$filename)) {
        fwrite(df, file = paste0(tempdir(), "/", get("filename", envir = envr), ".csv"))
    }
  
} # end HAC.simrep
