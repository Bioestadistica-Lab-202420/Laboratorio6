quality_ellipse <- function(X, alpha) {
  # Calculate mean and covariance
  result <- mean_cov_corr_matrices(X)
  M <- result$M
  VC <- result$VC
  
  # Eigen decomposition
  eig_result <- eigen(VC)
  eVect <- eig_result$vectors
  eVal <- eig_result$values
  eVal[c(1,2)] <- eVal[c(2,1)]
  eVal <- diag(eVal)
  eVect[ , c(1,2)]  <- eVect[ , c(2,1)]
  
  # Calculate axes lengths
  a <- sqrt(eVal[1, 1]) * sqrt(qchisq(1 - alpha, 2))
  b <- sqrt(eVal[2, 2]) * sqrt(qchisq(1 - alpha, 2))
  
  # Parametric equation
  t <- seq(0, 2 * pi, by = 0.01)
  nt <- length(t)
  x <- matrix(0, nrow = nt, ncol = 2)
  for (k in 1:nt) {
    x[k, ] <- (M + a * cos(t[k]) * eVect[, 1] + b * sin(t[k]) * eVect[, 2])
  }
  
  # Plot
  plot(x[, 1], x[, 2], type = "p", col = rgb(0.1, 0.3, 0.7, alpha = 0.5), pch = 19, cex = 0.4, xlab = "variable 1", ylab = "variable 2")
  points(M[1], M[2], col = "red", pch = 19, cex = 1.5)
  points(X[, 1], X[, 2], col = rgb(0.1, 0.3, 0.2, alpha = 0.5), pch = 19, cex = 0.4)
  legend("topright", legend = c("Mean", "Observations"), col = c("red", rgb(0.1, 0.3, 0.2, alpha = 0.5)), pch = 19, cex = 1.5)
  title(paste("Quality Ellipse at", (1 - alpha) * 100, "%"))
}
