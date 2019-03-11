\name{HACSim-package}
\alias{HACSim-package}
\alias{HACReal}
\alias{HACHypothetical}
\alias{HAC.simrep}
\docType{package}
\title{
\packageTitle{HACSim}
}
\description{
HACSim employs a novel method of iteratively generating species haplotype accumulation curves through extrapolation in order to assess sampling completeness based on the approach of Phillips et al. (2015).
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{HACSim}
\packageIndices{HACSim}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{HACSim}

Maintainer: \packageMaintainer{HACSim}
}
\references{
Phillips, J.D., Gwiazdowski, R.A., Ashlock, D. and Hanner, R. (2015). An exploration of sufficient sampling effort to describe intraspecific DNA barcode haplotype diversity: examples from the ray-finned fishes (Chordata: Actinopterygii). DNA Barcodes, 3: 66-73.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
## Simulate hypothetical species ##

N <- 100 # total number of sampled individuals
Hstar <- 10 # total number of haplotypes
probs <- rep(1/Hstar, Hstar) # equal haplotype frequency distribution

HACSObj <- HACHypothetical(N, Hstar, probs, filename = "output") # outputs a CSV file called "output.csv"
HAC.simrep(HACSObj) # runs a simulation


## Simulate real species ##

\dontrun{HACSObj <- HACReal(filename = "output") # outputs file called "output.csv"
HAC.simrep(HACSObj) # user prompted to select appropriate FASTA file of aligned/trimmed DNA sequences prior to running a simulation}
}
