# This function helps visualize the data in X. X is n by p matrix. 
# #The columns of X contains the observations of data points in the p variables
# The function is designed to explore the data in a 2D representation. 
# Therefore a scatter plot with 95# confidece ellipsoids is constructed to 
# visualize the data. A separate figure is presented for all combinatios of two 
# variables. The scatter plots shows association (R is shown on the plot) 
#
#INPUT X a n by p matrix. n observation on p variables
#
#OUTPUT Scatter plot for al combinations of variables n X two at a time. 

visualize_data <- function(X) {
  p <- ncol(X)
  
  source("./Funciones/mean_cov_corr_matrices.R")
  
  for (i in 1:p) {
    j <- i + 1
    while (j <= p){
      A <- cbind(X[, i], X[, j])
      result <- mean_cov_corr_matrices(A)
      M <- result$M
      VC <- result$VC
      R <- result$R[1, 2]
      
      e <- eigen(VC)
      eVect <- e$vectors
      eVal <- e$values
      eVal[c(1,2)] <- eVal[c(2,1)]
      eVal <- diag(eVal)
      
      eVect[ , c(1,2)]  <- eVect[ , c(2,1)]
      
      # Primera elipse
      m1 <- eVect[2, 1] / eVect[1, 1]
      theta1 <- atan(m1)
      x1 <- seq(M[1] - (sqrt(eVal[1, 1]) * sqrt(5.99)) * cos(theta1), M[1] + (sqrt(eVal[1, 1]) * sqrt(5.99)) * cos(theta1), length.out = 100)
      y1 <- m1 * x1 + M[2] - M[1] * m1
      
      # Segunda elipse
      m2 <- eVect[2, 2] / eVect[1, 2]
      theta2 <- atan(m2)
      x2 <- seq(M[1] - (sqrt(eVal[2, 2]) * sqrt(5.99)) * cos(theta2), M[1] + (sqrt(eVal[2, 2]) * sqrt(5.99)) * cos(theta2), length.out = 100)
      y2 <- m2 * x2 + M[2] - M[1] * m2
      
      # Gráfico
      plot(M[1], M[2], type = "p", pch = 19, col = "red", 
           main = paste("Scatter plot of X - R =", R),
           xlab = paste("X - ", i, " column"), 
           ylab = paste("X - ", j, " column"))
      
      points(A[, 1], A[, 2], pch = 19, col = "blue")
      lines(x1, y1, lty = 2, col = rgb(0, 0.4, 0.7, alpha = 0.8), lwd = 2)
      lines(x2, y2, lty = 2, col = rgb(0.4, 0, 0.7, alpha = 0.8), lwd = 2)
      
      # Añadir elipsoide
      center <- M
      a <- sqrt(eVal[1, 1]) * sqrt(5.99)
      b <- sqrt(eVal[2, 2]) * sqrt(5.99)
      t_vector <- seq(0, 2 * pi, 0.01)
      n <- length(t)
      ellipse <- matrix(0, nrow = n, ncol = 2)
      for (k in 1:n){
        ellipse[k,] <- t(center + a*cos(t_vector[k])*eVect[,1] + b*sin(t_vector[k])*eVect[,2])
      }
      #ellipse <- cbind(center[1] + a * cos(t) * eVect[, 1] + b * sin(t) * eVect[, 2], 
      #                 center[2] + a * sin(t) * eVect[, 1] - b * cos(t) * eVect[, 2])
      lines(ellipse[,1], ellipse[,2], col = rgb(0.1, 0.3, 0.2, alpha = 0.8), 
            pch = 19, lwd = 2)
      
      legend("topright", legend = c("Mean", "Data", "Axis lambda1", "Axis lambda2", "95% ellipsoid"), 
             col = c("red", "blue", rgb(0, 0.4, 0.7, alpha = 0.8), rgb(0.4, 0, 0.7, alpha = 0.8), rgb(0.1, 0.3, 0.2, alpha = 0.8)), 
             pch = 19, lty = c(1, 1, 2, 2, 1))
      
      plot_dimensions <- par("pin")
      
      #text(1.1 * min(A[, 1]), 0.75 * max(A[, 2]),  paste(VC, collapse = "\t"))
      
      for (i in 1:nrow(VC)) {
        for (j in 1:ncol(VC)) {
          text(1.1 * par("usr")[1] + (j-1)*12, 
               par("usr")[4] - 5 - (i-1)*8, 
               round(VC[i, j], digits = 4), cex = 1, font = 2)
        }
      }
      
      j <- j + 1
    }
  }
}
