onewayanova <- function(...) {
  np <- length(list(...)) # Number of populations
  # Organize observations in a vector
  y <- c() # Initialize y as an empty vector
  e <- c() # Initialize e (effects)
  xm <- numeric(np) # Means of each population
  n <- numeric(np) # Sizes of each population
  
  for (i in 1:np) {
    x <- unlist(list(...)[[i]])
    xm[i] <- mean(x) # Mean of the population
    n[i] <- length(x)
    y <- c(y, x) # Observations
    e <- c(e, rep(xm[i], n[i])) # Population means in a vector
  }
  nt <- sum(n)
  Vmean <- rep(mean(y), nt) # Vector with the overall mean
  Veffects <- e - Vmean # Vector with the differences between population mean and overall mean
  Vresiduos <- y - e # Vector with the residuals
  # Sum of squares
  SSeffects <- sum(Veffects^2)
  SSres <- sum(Vresiduos^2)
  # Degrees of freedom
  dfeff <- np - 1
  dfres <- nt - np
  # Statistic
  F <- (SSeffects/dfeff)/(SSres/dfres)
  # Critical F value with alpha = 0.05
  Fc <- qf(1 - 0.05, dfeff, dfres)
  # ANOVA Table
  AT <- data.frame(SS = c(SSeffects, SSres),
                   df = c(dfeff, dfres),
                   F = c(F, NA),
                   Fc = c(Fc, NA),
                   row.names = c('Efecto poblacional', 'Residuos'))
  return(AT)
}

# Example usage:
# AT <- onewayanova(X1, X2, X3)
# print(AT)
