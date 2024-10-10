confidence_region_two_populations <- function(X, Y) {
  # Calculate mean, covariance, and correlation matrices for X
  source("./funciones/mean_cov_corr_matrices.R")
  M1S1 <- mean_cov_corr_matrices(X)
  M1 <- M1S1$M
  S1 <- M1S1$VC
  
  # Calculate mean, covariance, and correlation matrices for Y
  M2S2 <- mean_cov_corr_matrices(Y)
  M2 <- M2S2$M
  S2 <- M2S2$VC
  
  # Get dimensions of X and Y
  n1 <- nrow(X)
  p <- ncol(X)
  n2 <- nrow(Y)
  
  alpha <- 0.05 # Confidence level
  
  # Calculate pooled covariance matrix
  Sp <- (n1 - 1) / (n1 + n2 - 2) * S1 + (n2 - 1) / (n1 + n2 - 2) * S2
  
  # Eigen decomposition
  eig_result <- eigen(Sp)
  eVect <- eig_result$vectors
  eVal <- eig_result$values
  
  eVal[c(1,2)] <- eVal[c(2,1)]
  eVal <- diag(eVal)
  eVect[ , c(1,2)]  <- eVect[ , c(2,1)]
  
  
  # Calculate critical T
  Tc <- qf(1 - alpha, p, n1 + n2 - p - 1) * p * (n1 + n2 - 2) / (n1 + n2 - p - 1)
  
  # Calculate center
  M <- M1 - M2
  
  # Add ellipsoid
  a <- sqrt(eVal[1,1]) * sqrt((1 / n1 + 1 / n2) * Tc) # Axis 1
  b <- sqrt(eVal[2,2]) * sqrt((1 / n1 + 1 / n2) * Tc) # Axis 2
  
  # Parametric equation
  t <- seq(0, 2 * pi, by = 0.01) # All angles
  nt <- length(t)
  x <- matrix(0, nrow = nt, ncol = 2)
  
  for (k in 1:nt) {
    x[k,] <- M + (a * cos(t[k]) * eVect[,1]) + (b * sin(t[k]) * eVect[,2])
  }
  
  # Plot
  plot(x[,1], x[,2], type = "p", col = "blue", pch = 20, cex = 0.4)
  points(M[1], M[2], col = "red", pch = 20, cex = 1.5)
  
  # Calculate simultaneous confidence intervals
  ul1 <- M1[1] - M2[1] + sqrt(Tc) * sqrt((1 / n1 + 1 / n2) * Sp[1, 1])
  ll1 <- M1[1] - M2[1] - sqrt(Tc) * sqrt((1 / n1 + 1 / n2) * Sp[1, 1])
  ul2 <- M1[2] - M2[2] + sqrt(Tc) * sqrt((1 / n1 + 1 / n2) * Sp[2, 2])
  ll2 <- M1[2] - M2[2] - sqrt(Tc) * sqrt((1 / n1 + 1 / n2) * Sp[2, 2])
  
  lines(seq(ll1, ul1, length.out = 100), rep(0, 100), lwd = 3, col = "yellow")
  lines(rep(0, 100), seq(ll2, ul2, length.out = 100), lwd = 3, col = "purple")
  
  title('Confidence interval for the mean')
}
