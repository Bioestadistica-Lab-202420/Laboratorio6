onewaymanova <- function(...) {
  # One way MANOVA construye la tabla de MANOVA al 0.05% para tantas poblaciones
  # como se pasen en los inputs de la función. 
  # Los datos de cada poblacion deben venir organizados en matrices de tamaño
  # p (numero de variables) por ng (numero de observaciones en cada
  # poblacion. El numero de observaciones puede ser distinto pero el número de
  # variables debe ser conservado en cada población. 
  # 
  # INPUT X1 Vector de tamaño p x n1 con observaciones de cada variable
  #          (filas) para la poblacion 1.
  # 
  #       Xg Matriz de tamaño p x n1 con observaciones de cada variable
  #          (filas) para la poblacion 1.
  # 
  # OUTPUT MT Manova Table 
  
  # Set the level of confidence
  alpha <- 0.05
  np <- length(list(...)) # Number of populations
  
  # Get information about the number of variables
  x1 <- list(...)[[1]] # Take the data from the first population
  p <- nrow(x1) # Calculate the number of variables (p) and the number of data in 1
  
  # Initialize variables for calculating means, S and total mean
  M <- matrix(0, nrow = p, ncol = np) # Initialize a matrix of population means (columns) by variable (rows)
  Y <- c() # Initialize an empty vector to save all the data together
  S <- vector("list", np) # Initialize a list to save the V-CV matrices of each population
  
  # Calculate mean and covariance
  source("./funciones/mean_cov_corr_matrices.R")
  
  # Get the mean and covariance matrix for each population
  for (i in 1:np) {
    x <- list(...)[[i]] # Take the data from the i-th population
    result <- mean_cov_corr_matrices(t(x))
    Mi <- result$M
    Si <- as.matrix(result$VC)
    S[[i]] <- Si
    M[, i] <- Mi # Save means
    Y <- rbind(Y, t(x)) # Save in Y all the concatenated data in a matrix n1+n2+... by p columns
  }
  
  # Calculate the total mean
  MT <- colMeans(Y)
  n <- nrow(Y)
  
  # Initialize matrices B and W
  B <- matrix(0, nrow = p, ncol = p)
  W <- matrix(0, nrow = p, ncol = p)
  N <- numeric(np)
  
  # Calculate sum of squares between and within
  for (i in 1:np) {
    x <- list(...)[[i]] # Take the data from the i-th population
    ni <- ncol(x)
    N[i] <- ni
    B <- B + ni * (M[, i] - t(t(MT))) %*% t(M[, i] - t(t(MT)))
    W <- W + (ni - 1) * S[[i]]
  }
  
  # Calculate the Wilks Lambda statistic
  Wlambda <- det(W) / det(B + W)
  
  MTab <- data.frame(matrix(ncol = 4, nrow = 0))
  
  # Depending on np (number of populations and the number of variables, the tests are different)
  if (p == 1 && np >= 2) {
    F <- (n - np) / (np - 1) * ((1 - Wlambda) / Wlambda)
    dfeff <- np - 1
    dfres <- n - np
    Fcritico <- qf(1 - alpha, dfeff, dfres)
    # MANOVA Table
    MTab <- data.frame(SS = c(det(W), det(B + W)),
                     df = c(dfeff, dfres),
                     F = c(F, NA),
                     Fc = c(Fcritico, NA),
                     row.names = c('Efecto poblacional', 'Residuos'))
  } else if (p == 2 && np >= 2) {
    F <- (n - np - 1) / (np - 1) * ((1 - sqrt(Wlambda)) / sqrt(Wlambda))
    dfeff <- 2 * (np - 1)
    dfres <- 2 * (n - np - 1)
    Fcritico <- qf(1 - alpha, dfeff, dfres)
    # MANOVA Table
    MTab <- data.frame(SS = c(det(W), det(B + W)),
                     df = c(dfeff, dfres),
                     F = c(F, NA),
                     Fc = c(Fcritico, NA),
                     row.names = c('Efecto poblacional', 'Residuos'))
  } else if (p >= 1 && np == 2) {
    F <- (n - p - 1) / (p - 1) * ((1 - Wlambda) / Wlambda)
    dfeff <- p - 1
    dfres <- n - p - 1
    Fcritico <- qf(1 - alpha, dfeff, dfres)
    # MANOVA Table
    MTab <- data.frame(SS = c(det(W), det(B + W)),
                     df = c(dfeff, dfres),
                     F = c(F, NA),
                     Fc = c(Fcritico, NA),
                     row.names = c('Efecto poblacional', 'Residuos'))
  } else if (p >= 1 && np == 3) {
    F <- (n - p - 2) / (p) * ((1 - sqrt(Wlambda)) / sqrt(Wlambda))
    dfeff <- 2 * p
    dfres <- 2 * (n - p - 2)
    Fcritico <- qf(1 - alpha, dfeff, dfres)
    # MANOVA Table
    MTab <- data.frame(SS = c(det(W), det(B + W)),
                     df = c(dfeff, dfres),
                     F = c(F, NA),
                     Fc = c(Fcritico, NA),
                     row.names = c('Efecto poblacional', 'Residuos'))
  } else {
    Chi <- -(n - 1 - (p + np) / 2) * log(Wlambda)
    dfeff <- p * (np - 1)
    Chicritico <- qchisq(1 - alpha, dfeff)
    # MANOVA Table
    MTab <- data.frame(SS = c(det(W), det(B + W)),
                     df = c(dfeff, NA),
                     Chi = c(Chi, NA),
                     Chic = c(Chicritico, NA),
                     row.names = c('Efecto poblacional', 'Residuos'))
  }
  return(list(MT = MTab, M = M, N = N, W = W))
}

# Example usage:
# MT <- onewaymanova(X1, X2, X3)
# print(MT)
