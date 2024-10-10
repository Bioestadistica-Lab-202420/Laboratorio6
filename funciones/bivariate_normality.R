bivariate_normality <- function(X) {
  # detecta normalidad para una y 2+ variables
  # chequear si X tiene una sola columna o varias
  if ("data.frame" %in% class(X)){
    X <- as.matrix(X)
  }
  
  nc <- ncol(X)
  nr <- nrow(X)
  if (nc == 1) {
    # QQ plot
    xsorted <- sort(X)
    j <- 1:nr
    samplequartil <- (j - 0.5) / nr
    q <- qnorm(samplequartil)
    T <- data.frame(xsorted = xsorted, samplequartil = samplequartil, q = q)
    nvec <- NULL
  } else {
    # Bivariate
    source("./Funciones/mean_cov_corr_matrices.R")
    result <- mean_cov_corr_matrices(X)
    M <- result$M
    print(dim(M))
    VC <- result$VC
    print(VC)
    SI <- solve(VC)
    print(dim(SI))
    GD <- numeric(nr)
    for (i in 1:nr) {
      GD[i] <- t((X[i,]) - M) %*% SI %*% ((X[i,]) - M)
    }
    vc <- qchisq(0.5, nc)
    nvec <- (nr - sum(GD > vc)) / nr
    GDsorted <- sort(GD)
    ChiSquare <- numeric(nr)
    for (i in 1:nr) {
      ChiSquare[i] <- qchisq((i - 0.5) / nr, nc)
    }
    comparison <- GDsorted > vc
    Observations <- 1:nr
    T <- data.frame(Observations = Observations, GDsorted = GDsorted, ChiSquare = ChiSquare, comparison = comparison)
  }
  return(list(T = T, nvec = nvec))
}