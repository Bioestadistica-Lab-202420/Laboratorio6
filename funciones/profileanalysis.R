profileanalysis <- function(X1, X2, q) {
  # Define test reliability
  alpha <- 0.05
  # Numbers needed for the tests
  g <- 2
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)
  # Means and Var-Cov matrices for the groups
  M1 <- colMeans(X1)
  S1 <- cov(X1)
  M2 <- colMeans(X2)
  S2 <- cov(X2)
  # Test for parallelism
  C <- matrix(0, nrow = p - 1, ncol = p)
  for (i in 1:(p - 1)) {
    for (j in 1:p) {
      if (i == j) {
        C[i, j] <- -1
      } else if (j == i + 1) {
        C[i, j] <- 1
      }
    }
  }
  N <- n1 + n2
  W <- (n1 - 1) * S1 + (n2 - 1) * S2
  Spooled <- (1 / (N - g)) * W
  T <- t(M1 - M2) %*% t(C) %*% solve((1 / n1 + 1 / n2) * C %*% Spooled %*% t(C)) %*% C %*% (M1 - M2)
  Tcritico <- (n1 + n2 - 2) * (p - 1) / (n1 + n2 - p) * qf(1 - alpha, p - 1, n1 + n2 - p)
  if (T < Tcritico) {
    cat('Los grupos tienen un perfil paralelo. Continua el analisis de coincidencia\n')
    O <- matrix(1, nrow = p, ncol = 1)
    Tcon <- t(O) %*% (M1 - M2) %*% solve((1 / n1 + 1 / n2) * t(O) %*% Spooled %*% O) %*% t(O) %*% (M1 - M2)
    Tconcritico <- qf(1 - alpha, 1, n1 + n2 - 2)
    if (Tcon < Tconcritico) {
      cat('Los grupos tienen un perfil coincidente. No hay evidencia de efecto del tratamiento. Continua el analisis de horzontalidad del perfil\n')
      MP <- n1 / (n1 + n2) * M1 + n2 / (n1 + n2) * M2
      Thor <- (n1 + n2) * t(MP) %*% t(C) %*% solve(C %*% Spooled %*% t(C)) %*% C %*% MP
      Thorcritico <- (n1 + n2 - 1) * (p - 1) / (n1 + n2 - p + 1) * qf(1 - alpha, p - 1, n1 + n2 - p + 1)
      if (Thor < Thorcritico) {
        cat('Los grupos tienen un perfil horizontal. Termina el analisis de perfil\n')
      } else {
        cat('Los grupos no tienen un perfil horizontal. Se procede a evaluar la hipotesis del polinomio del grado suministrado\n')
        B <- matrix(0, nrow = p, ncol = q + 1)
        for (i in 0:(p - 1)) {
          for (j in 0:q) {
            B[i + 1, j + 1] <- i^j
          }
        }
        BetaC <- solve(t(B) %*% solve(Spooled) %*% B) %*% t(B) %*% solve(Spooled) %*% M1
        BetaT <- solve(t(B) %*% solve(Spooled) %*% B) %*% t(B) %*% solve(Spooled) %*% M2
        WC <- matrix(0, nrow = p, ncol = p)
        for (j in 1:n1) {
          WC <- WC + (as.matrix(X1[j,]) - B %*% BetaC) %*% t(as.matrix(X1[j,]) - B %*% BetaC)
        }
        WT <- matrix(0, nrow = p, ncol = p)
        for (j in 1:n2) {
          WT <- WT + (as.matrix(X2[j,]) - B %*% BetaT) %*% t(as.matrix(X2[j,]) - B %*% BetaT)
        }
        Wq <- WC + WT
        WLambda <- det(W) / det(Wq)
        Chi <- -(N - 1 / 2 * (p - q + g)) * log(WLambda)
        Chicritico <- qchisq(1 - alpha, (p - q - 1) * g)
        if (Chi < Chicritico) {
          cat('El polinomio ajusta de manera significativa\n')
        } else {
          cat('El polinomio no ajusta\n')
        }
      }
    } else {
      cat('Los grupos tienen un perfil no coincidente. Hay evidencia de efectos del tratamiento\n')
    }
  } else {
    cat('Los grupos tienen un perfil no paralelo. No se continua con el analisis de concordancia\n')
  }
  # Graphical output
  plot(M1, type='o', pch=19, cex=1, col='blue', xlab='mediciones', ylab='Respuesta',
       ylim = c(min(c(M1,M2)), max(c(M1,M2))))
  points(M2, type='o', pch=19, cex=1, col='orange')
  legend('topright', legend=c('X1', 'X2'), pch=19, col=c('blue', 'orange'))
}