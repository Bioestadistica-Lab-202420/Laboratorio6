FactorAnalysis <- function(LL, LH, HL, HH) {
  # Define parameters
  g <- 2  # Number of factors
  b <- 2  # Number of levels
  n <- nrow(LL)  # Number of observations
  p <- ncol(LL)  # Number of variables
  
  # Determine means and covariance matrices
  source("./funciones/mean_cov_corr_matrices.R")
  
  result <- mean_cov_corr_matrices(LL)
  MLL <- result$M
  SLL <- as.matrix(result$VC)
  
  result <- mean_cov_corr_matrices(LH)
  MLH <- result$M
  SLH <- as.matrix(result$VC)
  
  result <- mean_cov_corr_matrices(HL)
  MHL <- result$M
  SHL <- as.matrix(result$VC)
  
  result <- mean_cov_corr_matrices(HH)
  MHH <- result$M
  SHH <- as.matrix(result$VC)
  
  # Means for each level of each factor
  MF1L <- colMeans(rbind(LL, LH))
  MF1H <- colMeans(rbind(HL, HH))
  MF2L <- colMeans(rbind(LL, HL))
  MF2H <- colMeans(rbind(LH, HH))
  
  
  # Overall mean
  M <- colMeans(rbind(LL, LH, HL, HH))
  
  # Sum of squares
  SSF1 <- b * n * tcrossprod(MF1L - M) + b * n * tcrossprod(MF1H - M)
  SSF2 <- g * n * tcrossprod(MF2L - M) + g * n * tcrossprod(MF2H - M)
  SSF1F2 <- n * (tcrossprod(MLL - MF1L - MF2L + M) + 
                   tcrossprod(MLH - MF1L - MF2H + M) + 
                   tcrossprod(MHL - MF1H - MF2L + M) + 
                   tcrossprod(MHH - MF1H - MF2H + M))
  SSRes <- (n - 1) * (SLL + SLH + SHL + SHH)
  
  # Calculate Wilks' lambda
  Wlambda <- det(SSRes) / det(SSRes + SSF1F2)
  
  # Calculate F statistic
  F <- ((1 - Wlambda) / Wlambda) * ((g * b * (n - 1) - p + 1) / 2) / ((abs((g - 1) * (b - 1) - p) + 1) / 2)
  
  # Calculate critical F value
  Fc <- qf(1 - 0.05, abs((g - 1) * (b - 1) - p) + 1, g * b * (n - 1) - p + 1)
  
  # Calculate F for each factor
  WlambdaF1 <- det(SSRes) / det(SSRes + SSF1)
  WlambdaF2 <- det(SSRes) / det(SSRes + SSF2)
  F1 <- ((1 - WlambdaF1) / WlambdaF1) * ((g * b * (n - 1) - p + 1) / 2) / ((abs((g - 1) * (b - 1) - p) + 1) / 2)
  F2 <- ((1 - WlambdaF2) / WlambdaF2) * ((g * b * (n - 1) - p + 1) / 2) / ((abs((g - 1) * (b - 1) - p) + 1) / 2)
  
  # Create MANOVA table
  MTF <- data.frame(SS = c(det(SSRes + SSF1F2), det(SSRes + SSF1), det(SSRes + SSF2), det(SSRes)),
                    df = c((g - 1) * (b - 1), g - 1, b - 1, g * b * (n - 1)),
                    F = c(F, F1, F2, NA),
                    Fc = c(Fc, NA, NA, NA))
  rownames(MTF) <- c('Interaccion', 'Factor 1', 'Factor 2', 'Residuos')
  return(MTF)
}
