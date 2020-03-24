\name{sim.seqs}
\alias{sim.seqs}

\title{Simulate DNA sequences according to DNA substitution models}

\description{}

\usage{sim.seqs(num.seqs, num.haps, length.seqs, nucl.freqs, count.haps, subst.model = c("JC69", "K80", "F81", "HKY85"), mu.rate, transi.rate, transv.rate)}

\arguments{\item{num.seqs}{Number of simulated DNA sequences}
  \item{num.haps}{Number of simulated unique species' haplotypes}
    \item{length.seqs}{Length of DNA sequences}
    \item{nucl.freqs}{Nucleotide frequency distribution vector of A, C, G, and T respectively}
    \item{count.haps}{Haplotype frequency distribution vector}
    \item{subst.model}{Model of DNA substitution}
    \item{mu.rate}{Overall nucleotide mutation rate}
    \item{transi.rate}{Nucleotide transition rate}
    \item{transv.rate}{Nucleotide transversion rate}
  }
    
    \value{A FASTA file of DNA sequences}
    
\note{\code{num.seqs} must be greater than or equal to \code{num.haps}.

  Both \code{num.seqs} and \code{num.haps} must be greater than 1.

  \code{nucl.freqs} must have a length of four and its elements must sum to 1.
  
  \code{count.haps} must have a length of \code{num.haps} and its elements must sum to \code{num.seqs}. 
  
  \code{subst.model} must be one of "JC69" (Jukes Cantor corrected p-distance), "K80" (Kimura-2-Parameter (K2P), "F81" (Felenstein) or "HKY85"(Hasegawa-Kishino-Yano)
  
  \code{mu.rate} must be specified for both "JC69" and "F81" models
  
  \code{transi.rate} and \code{transv.rate} must be specified for both "K80" and "HKY85" models
  
  All elements \code{nucl.freqs} must be equal to 0.25 when \code{subst.model} is either "JC69" or "K80" 
  
  All elements \code{nucl.freqs} must differ from 0.25 when \code{subst.model} is either "F81" or "HKY85"
  
  }

\examples{

}