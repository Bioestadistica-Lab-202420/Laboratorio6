Tsquare_chart <- function(X, alpha) {
  # Calculate mean and covariance
  result <- mean_cov_corr_matrices(X)
  M <- result$M
  VC <- result$VC
  SI <- solve(VC)
  
  # Calculate T2 statistic
  n <- nrow(X)
  p <- ncol(X)
  GD <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    GD[i] <- t(X[i, ] - M) %*% SI %*% (X[i, ] - M)
  }
  
  # Calculate critical value
  vc <- qchisq(1 - alpha, p)
  
  # Plot
  plot(1:n, rep(vc, n), type = "l", lty = 2, col = "black", xlab = "observation", 
       ylab = "T2 statistic", 
       main = paste("T Quality Chart", (1 - alpha) * 100, "%"),
       xlim = c(0,n), ylim = c(min(GD),max(max(GD),vc)))
  points(1:n, GD, col = rgb(0.1, 0.3, 0.2, alpha = 0.5), pch = 19, cex = 0.4)
}