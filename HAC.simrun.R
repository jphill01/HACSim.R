##### Set working directory #####

setwd("/Users/jarrettphillips/Desktop/HAC simulation")
# setwd("/Users/jarrettphillips/Desktop/HACSim-RShiny-App-master")

##### Clear memory #####

remove(list = ls())

##### Load packages #####

library(HACSim)
?HACSim


##### SCENARIOS #####

# N = 10, 50, 100
# Hstar = 10, 20, 25

## Example probs for Hstar = 10 - you can come up with others ##

# probs = rep(1/Hstar, Hstar) # equal haplotype frequency
# probs = c(0.90, rep(0.10/9, 9)) # 1 dominant haplotype
# probs = c(0.45, 0.45, rep(0.10/8, 8)) # 2 dominant haplotypes
# probs = c(0.30, 0.30, 0.30, rep(0.10/7, 7)) # 3 dominant haplotypes

## Haplotype Recovery ## 

# p = 0.90
# p = 0.95 - default
# p = 0.99
# p = 1 - long runtime


### SET PARAMETERS ###

## Set parameters for hypothetical species ##

N <- 27 # total number of sampled individuals
Hstar <- 11 # total number of haplotypes
# probs <- c(rep(0.3, 3), rep(0.1/42, 42)) # must sum to 1
# probs <- rep(1/Hstar, Hstar) # equal haplotype frequency
probs <- c(7, 5, 4, rep(2, 3), rep(1, 5)) / N
# probs <- c(53, 13, 9, 1) / N
# probs <- c(96, 2) / N
# probs <- c(215, rep(3, 2), rep(2, 2), rep(1, 10)) / N
# probs <- c(0.3, 0.3, 0.3, rep(0.10/5, 5))
# probs <- c(0.45, 0.45, rep(0.10/8, 8))
# probs <- c(0.9, 0.025, 0.025, 0.025, 0.025)
N <- envr$N
Hstar <- envr$Hstar
probs <- envr$probs


### Run HAC Simulations ###

## Simulate hypothetical species ##
## If not set, prop defaults to 0.10
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 100, p = 0.95, conf.level = 0.95, num.iters = NULL, filename = NULL)

## Simulate hypothetical species - subsampling ##
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 10000, p = 0.95, subsample = TRUE, prop = 0.25, ci.type = "asymptotic", conf.level = 0.95, num.iters = NULL, filename = "output")

## Simulate hypothetical species and all parameters changed - subsampling ##
HACSObj <- HACHypothetical(N = N, Hstar = Hstar, probs = probs, perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, ci.type = "asymptotic", conf.level = 0.95, filename = "output")

## Simulate real species ##
HACSObj <- HACReal(p = 0.95, perms = 50, ci.type = "asymptotic", conf.level = 0.95, num.iters = 1, filename = "output")

## Simulate real species - subsampling ##
HACSObj <- HACReal(p = 0.95, subsample = TRUE, prop = 0.25, ci.type = "asymptotic", conf.level = 0.95, num.iters = NULL, filename = "output")

## Simulate real species and all parameters changed - subsampling ##
HACSObj <- HACReal(perms = 10000, p = 0.90, subsample = TRUE, prop = 0.15, ci.type = "asymptotic", conf.level = 0.99, filename = "output")

## Set parameters to simulate DNA sequences ##

num.seqs <- 27 # number of DNA sequences
num.haps <- 11 # number of haplotypes
length.seqs <- 658 # length of DNA sequences
count.haps <- c(7, 5, 4, rep(2, 3), rep(1, 5)) # haplotype frequency distribution
nucl.freqs <- c(0.3335588, 0.1191161, 0.1052483, 0.4420768 ) # nucleotide frequency distribution
subst.model <- "HKY85" # desired nucleotide substitution model
codon.tbl <- "vertebrate mitochondrial"
mu.rate <- NULL # mutation rate
transi.rate <- 1e-2 # transition rate
transv.rate <- transi.rate / 2 # transversion rate

sim.seqs(num.seqs = num.seqs, num.haps = num.haps, length.seqs = length.seqs, nucl.freqs = nucl.freqs, count.haps = count.haps, codon.tbl = codon.tbl, subst.model = subst.model, mu.rate = mu.rate, transi.rate = transi.rate, transv.rate = transv.rate)

HACSObj <- HACReal(p = 0.95, perms = 10000, ci.type = "asymptotic", conf.level = 0.95, num.iters = NULL, filename = "output")

### Run HAC simulations ###

# set.seed(0673227) # set random seed for reproducibility, if desired

# The below loop works correctly

reps <- 100
res1 <- vector(reps, mode = "numeric") # preallocate vector
# res2 <- vector(reps, mode = "numeric") # preallocate vector
ptm <- proc.time()
for (i in 1:reps) {
  HAC.simrep(HACSObj) # simulate
  res1[i] <- envr$Nstar - envr$X # store values of N*
  # res2[i] <- envr$R # store values of haplotype recovery 
  }
proc.time() - ptm


plot(sort(table(res1)), xlab = "N*", ylab = "Frequency")
mean(res1) 

par(mfrow = c(2, 2))

hist(res1)
abline(v = mean(res1), lwd = 2, col = "red")
plot(res1, type = "l")
abline(h = mean(res1), lwd = 2, col = "red")

summary(res1)
quantile(res1, c(0.025, 0.975))
tail(sort(table(res1)), n = 1) # most frequent N* value

summary(res2)
quantile(res2, c(0.025, 0.975))
tail(sort(table(signif(res2, digits = 2))), n = 1) # most frequent R value


# Iteration plots - can only be plotted for a single run

df_out <- as.matrix(envr$df_out)
df_out <- matrix(df_out, ncol = ncol(df_out), dimnames = NULL)

plot(1:envr$iters, df_out[, 5] / df_out[,1], type = "l", xlab = "Number of iterations", ylab = "Step size")
plot(1:envr$iters, (envr$Nstar - df_out[,5]) / (envr$Hstar - df_out[,1] + 1e6), type = "l", xlab = "Number of iterations", ylab = "Step size")


par(mfrow = c(2, 3))

plot(1:envr$iters, df_out[, 1], type = "l", xlab = "Number of iterations", ylab = "Mean number of haplotypes sampled") # number of iterations vs. mean number of haplotypes sampled
plot(1:envr$iters, df_out[, 2], type = "l", xlab = "Number of iterations", ylab = "Mean number of haplotypes not sampled") # number of iterations vs. mean number of haplotypes not sampled
plot(1:envr$iters, df_out[, 3], type = "l", xlab = "Number of iterations", ylab = "Proportion of haplotypes sampled") # number of iterations vs. proportion of haplotypes sampled
plot(1:envr$iters, df_out[, 4], type = "l", xlab = "Number of iterations", ylab = "Proportion of haplotypes not sampled") # number of iterations vs. proportion of haplotypes not sampled
plot(1:envr$iters, df_out[, 5], type = "l",  xlab = "Number of iterations", ylab = "Mean value of N*") # number of iterations vs. mean value of N*
plot(1:envr$iters, df_out[, 6], type = "l", xlab = "Number of iterations", ylab = "Number of additional specimens required to be sampled") # number of iterations vs. number of additional specimens required to be sampled


