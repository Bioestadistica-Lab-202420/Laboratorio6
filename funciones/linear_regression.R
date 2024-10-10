linear_regression <- function(Y, X) {
  library(pracma)
  alpha <- 0.05
  n <- length(Y)
  r <- ncol(X)
  
  # Calculate mean of Y
  Ym <- mean(Y)
  
  # Design matrix
  Z <- as.matrix(cbind(rep(1, n), X))
  
  # Compute (Z'Z)^-1
  W <- pinv(t(Z) %*% Z)
  
  # Compute Betas
  Beta <- W %*% t(Z) %*% Y
  
  # Compute Y estimated
  Yh <- Z %*% Beta
  
  # Compute error (Residuals)
  E <- Y - Yh
  
  # Residuals sum of squares
  SSerror <- sum(E^2)
  
  # Compute s^2 (variance)
  df <- n - r - 1
  s2 <- SSerror / df
  
  # Coefficient of determination
  SSregresion <- sum((Yh - Ym)^2)
  SStotal <- sum((Y - Ym)^2)
  R <- SSregresion / SStotal
  
  # Hat matrix
  H <- Z %*% W %*% t(Z)
  
  # Studentized errors
  Eest <- E / sqrt(s2 * (1 - diag(H)))
  
  if (!all(X == 0 | X == 1)) {
  #if (all(X != 0) && all(X != 1)) {
    # Produce variance-covariance matrix for elements in W
    VAR <- s2 * W
    
    # Confidence intervals for the betas
    IC <- matrix(0, nrow = r + 1, ncol = 2)
    betas <- list()
    for (i in 1:(r + 1)) {
      IC[i, 1] <- Beta[i] - sqrt(VAR[i, i]) * sqrt((r + 1) * qf(1 - alpha, r + 1, df))
      IC[i, 2] <- Beta[i] + sqrt(VAR[i, i]) * sqrt((r + 1) * qf(1 - alpha, r + 1, df))
      betas[[i]] <- paste("Beta", i)
    }
    lower_bound <- IC[, 1]
    upper_bound <- IC[, 2]
    Beta_estimate <- Beta
    BT <- data.frame(Betas = unlist(betas), Beta_estimate, lower_bound, upper_bound)
    
    # Evaluation of the variables in the model
    if (r > 1) {
      # Test for each of the independent variables
      Fc <- qf(1 - alpha, r - 1, df)
      F <- numeric(r)
      variables <- c()
      suggestion <- c()
      for (i in 1:r) {
        Zi <- Z
        Zi <- Zi[, -(i + 1)] 
        Wi <- solve(t(Zi) %*% Zi)
        Betai <- Wi %*% t(Zi) %*% Y
        Yhi <- Zi %*% Betai
        Ei <- Y - Yhi
        SSerrori <- sum(Ei^2)
        F[i] <- ((SSerrori - SSerror) / (r - 1)) / s2
        variables <- append(variables,paste("Variable", i))
        if (F[i] > Fc) {
          suggestion <- append(suggestion,"include")
        } else {
          suggestion <- append(suggestion,"exclude")
        }
      }
      F_critico <- rep(Fc, r)
      VT <- data.frame(variables, F_test = F, F_critico, suggestion)
    } else {
      VT <- c()
      cat("Only one variable, no hold out procedure available for testing significance of independent variables.\n")
    }
    
    # Output plot
    if (r == 1) {
      plot(X, Y, type = "p", col = "red", pch = 19, cex = 1.5, xlab = "X", ylab = "Y")
      lines(X, Yh, type = "l", col = "blue", lwd = 1)
    }
    AT <- c()
  } else {
    SourceVariation <- c("Model", "Error", "Total")
    SumSquare <- c(SSregresion, SSerror, SStotal)
    DegFree <- c(r + 1, n - r, n - 1)
    MeanSquare <- c(SSregresion / (r + 1), SSerror / (n - r), 0)
    FValue <- c(MeanSquare[1] / MeanSquare[2], 0, 0)
    Probablity <- c(pf(FValue[1], DegFree[1], DegFree[2]), 0, 0)
    AT <- data.frame(SourceVariation, SumSquare, DegFree, MeanSquare, FValue, Probablity)
    BT <- c()
    VT <- c()
  }
  
  return(list(BT = BT, VT = VT, AT = AT, Yh = Yh, Eest = Eest))
}