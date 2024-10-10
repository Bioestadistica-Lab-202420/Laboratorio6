power_transformation <- function(X) {
  # ensayar varios valores de lambda 
  l <- seq(-2, 2, by = 0.01)
  nl <- length(l)
  nc <- ncol(X)
  nr <- nrow(X)
  F_m <- matrix(0, nrow = nl, ncol = nc)
  lambda <- numeric(nc)
  fmax <- numeric(nc)
  for (p in 1:nc) { #para cada columna de X
    for (i in 1:nl) { #para cada valor de lambda
      x <- X[, p] #seleccione la variable
      if (l[i] != 0) {
        xl <- ((x^l[i]) - 1) / l[i] #calcule la potencia (x^l-1)/l
      } else {
        xl <- log(x) #calcule cuando l es 0 
      }
      xlm <- mean(xl) #saque la media 
      sl <- (xl - xlm)^2 #saque la varianza
      F_m[i, p] <- -nr/2 * log(1/nr * sum(sl)) + (l[i] - 1) * sum(log(x)) #calcula la funcion de lambda
    }
    ix <- which.max(F_m[, p])
    lambda[p] <- l[ix]
    fmax[p] <- max(F_m[, p])
  }
  return(list(l = l, F_m = F_m, lambda = lambda, fmax = fmax))
}