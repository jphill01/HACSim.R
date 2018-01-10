### Data Type Conversion Functions ###

array2matrix <- function(pop) {
  d <- dim(pop)
  ind <- seq_along(d)
  pop <- aperm(pop, c(ind[2], ind[-2]))
  dim(pop) <- c(d[2], prod(d[-2]))
  pop <- t(pop)
  attr(pop, "origdim") <- d
  pop
}

matrix2array <- function(pop, d = attr(pop, "origdim")) {
  ind <- seq_along(d)
  pop <- t(pop)
  dim(pop) <- c(d[2], d[-2])
  aperm(pop, c(ind[2], ind[-2]))
}