## Master class, both functions return an object of this type
HACClass <- function(input.seqs = NULL,
                     sim.seqs = NULL,
                     num.seqs = NULL,
                     length.seqs = NULL,
                     nucl.freq = NULL,
                     subst.model = "JC69",
                     mu.rate = NULL,
                     transi.rate = NULL,
                     transv.rate = NULL,
                     subset.seqs = NULL,
                     prop.seqs = NULL,
                     prop.haps = NULL,
                     subset.haps = NULL,
                     N = NA,
                     Hstar = NA,
                     probs = NA,
                     p = NA,
                     perms = NA,
                     filename = NULL) {
    HACObject <- list(input.seqs = input.seqs,
                      sim.seqs = sim.seqs,
                      num.seqs = num.seqs,
                      length.seqs = length.seqs,
                      nucl.freq = nucl.freq,
                      subst.model = subst.model,
                      mu.rate = mu.rate,
                      transi.rate = transi.rate,
                      transv.rate = transv.rate,
                      subset.seqs = subset.seqs,
                      prop.seqs = prop.seqs,
                      prop.haps = prop.haps,
                      subset.haps = subset.haps,
                      N = N,
                      Hstar = Hstar,
                      probs = probs,
                      p = p,
                      perms = perms,
                      filename = filename
    )
    
    ## Set the name for the class
    class(HACObject) <- append(class(HACObject), "HAC")
    return(HACObject)
}

## This creates an evolution scenario or uses raw data from a file
HAC <- function(type,
                perms = 10000,
                p = 0.95,
                subsample = FALSE,
                prop = 0.1,
                filename = NULL) {
    ## Type has to be set
    if (type == "evolution") {
        # input.seqs <- FALSE # analyze DNA sequence file? 
        # subset.haps <- NULL # subset haplotypes? 
        # prop.haps <- NULL # proportion of haplotypes to subsample 
        # subset.seqs <- FALSE # subset DNA sequences? 
        # prop.seqs <- NULL # proportion of DNA sequences to subsample 
        # sim.seqs <- TRUE # simulate DNA sequences? 
        # subst.model <- "JC69" # nucleotide substitution model
        # prop.seqs <- prop # proportion of DNA sequences to subsample

        objectHAC <- HACClass(input.seqs = FALSE, 
                              subset.seqs = FALSE, 
                              sim.seqs = TRUE,
                              num.seqs = num.seqs,
                              length.seqs = length.seqs,
                              nucl.freq = NULL,
                              subst.model = "JC69",
                              mu.rate = mu.rate,
                              transi.rate = transi.rate,
                              transv.rate = transv.rate,
                              perms = perms, 
                              p = p,
                              filename = filename)
    } else if (type == "real") {
        # input.seqs <- TRUE # analyze DNA sequence file? 
        # sim.seqs <- FALSE # simulate DNA sequences? 
        # subset.haps <- NULL # subset haplotypes?  
        # prop.haps <- NULL # proportion of haplotypes to subsample 
        # subset.seqs <- TRUE # subset DNA sequences? 
        
        if (subsample == TRUE) {
            subset.seqs <- TRUE
            prop.seqs <- prop # proportion of DNA sequences to subsample
        } else {
            subset.seqs <- FALSE
            prop.seqs <- NA # proportion of haplotypes to subsample 
        }
        
        objectHAC <- HACClass(input.seqs = TRUE,
                              subset.seqs = subset.seqs,
                              sim.seqs = FALSE,
                              prop.seqs = prop.seqs,
                              perms = perms,
                              p = p,
                              filename = filename)
    }
    
    return(objectHAC)
}

## For simulating custom data distributions
HACHypothetical <- function(N,
                            Hstar,
                            probs,
                            perms = 10000,
                            p = 0.95,
                            subsample = FALSE,
                            prop = 0.1,
                            filename = NULL) {
    if (missing(N)) {
      stop("Please provide a value for N")
    }
    if (missing(Hstar)) {
        stop("Please provide a value for Hstar")
    }
    if (missing(probs)) {
        stop("Please provide a value for probs")
    }
    
    # Representation of what the parameters will be
  
    # input.seqs <- FALSE # subset DNA sequences? 
    # sim.seqs <- FALSE # simulate DNA sequences? 
    # subset.seqs <- FALSE # subset DNA sequences? 
    # prop.seqs <- NULL # proportion of DNA sequences to subsample 
    # subset.haps <- NULL # subset haplotypes? 
    
    if (subsample == TRUE) {
        prop.haps <- prop # proportion of haplotypes to subsample 
        subset.haps <- sort(sample(Hstar, size = ceiling(prop.haps * Hstar), replace = FALSE))
    } else  {
        prop.haps <- NULL # proportion of haplotypes to subsample
        subset.haps <- NULL
    }
  objectHAC <- HACClass(N = N,
                        Hstar = Hstar,
                        probs = probs,
                        input.seqs = FALSE,
                        subset.seqs = FALSE,
                        sim.seqs = FALSE,
                        prop.haps = prop.haps,
                        subset.haps = subset.haps,
                        perms = perms,
                        p = p,
                        filename = filename)
    return(objectHAC)
}
