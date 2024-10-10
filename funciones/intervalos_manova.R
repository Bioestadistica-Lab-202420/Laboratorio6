intervalos_manova <- function(M, N, W, Contraste) {
  # All necessary inputs for this function are produced by onewaymanova.
  # INPUT:
  # M: vector of population means where rows represent variables and columns represent populations (p rows, g columns)
  # N: vector of sample sizes for each population (1 row, g columns)
  # W: matrix of SS within
  # Contraste: a vector indicating the populations to compare. The algorithm returns the interval for all variables of that population combination. Example [1 2] compares population 1 with 2.
  
  alpha <- 0.05
  p1 <- Contraste[1]
  p2 <- Contraste[2]
  n1 <- N[p1]
  n2 <- N[p2]
  np <- length(N) # Number of populations
  p <- nrow(M) # Number of variables
  
  # Interval
  LI <- M[, p1] - M[, p2] - qt(1 - (alpha / (p * np * (np - 1))), sum(N) - np) * sqrt(diag(W) / (sum(N) - np) * (1 / n1 + 1 / n2))
  UI <- M[, p1] - M[, p2] + qt(1 - (alpha / (p * np * (np - 1))), sum(N) - np) * sqrt(diag(W) / (sum(N) - np) * (1 / n1 + 1 / n2))
  
  return(list(LI = LI, UI = UI))
}

# Example usage:
# intervalos_manova(M, N, W, c(1, 2))
