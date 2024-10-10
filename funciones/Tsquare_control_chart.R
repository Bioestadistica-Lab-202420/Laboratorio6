Tsquare_control_chart <- function(X, alpha) {
  # Calculate mean and covariance
  result <- mean_cov_corr_matrices(X)
  M <- result$M
  VC <- result$VC
  SI <- solve(VC)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Calculate T2 statistic
  GD <- rep(0, n)
  for (i in 1:n) {
    GD[i] <- (n / (n + 1)) * t(X[i, ] - M) %*% SI %*% (X[i, ] - M)
  }
  
  # Calculate critical value
  vc <- qf(1 - alpha, p, n - p) * p * (n - 1) / (n - p)
  
  # Plot
  plot(1:n, rep(vc, n), type = "l", lty = 2, col = "black", xlab = "observation",
       ylab = "T2 statistic",
       main = paste("T Control Chart", (1 - alpha) * 100, "%"),
       xlim = c(0,n), ylim = c(min(GD),max(max(GD),vc)))
  points(1:n, GD, col = rgb(0.1, 0.3, 0.2, alpha = 0.5), pch = 19, cex = 0.4)
}