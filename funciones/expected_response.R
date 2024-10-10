expected_response <- function(Y, X, Xo, type) {
  # la funcion exepcted response acepta un vector (univariado) de respuesta Y y
  # una matriz X de n por r donde n es el numero de observaciones y r es el
  # numero de variables independientes presentadas para el modelo de regresion lineal. 
  # Adicionamente recibe un vector de observaciones (hipoteticas (forecast) o
  # medidas (valor esperado)) y produce un valor de respuestas (expected
  # response) 
  #
  # INPUT Y Vector de n por 1 con las variables de respuesta
  #      X Matriz de n por r con las variables independientes. 
  #      Xo Vector con observaciones (1*r)
  #      type 1 = expected values, 2 = forecast
  #
  # OUPUT  Una tabla con el valor esperado o forecast de Y y sus intervalos de
  #        confianza 
  
  alpha <- 0.05
  n <- nrow(X)
  r <- ncol(X)
  # design matrix
  Z <- as.matrix(cbind(1, X))
  
  # compute (Z'Z)^-1
  W <- solve(t(Z) %*% Z)
  # compute Bs
  Beta <- W %*% t(Z) %*% Y
  # compute Y estimated
  Yh <- Z %*% Beta
  # compute error (Residuals)
  E <- Y - Yh
  # residuals sum of squares
  SSerror <- t(E) %*% E
  # cumpute s^2 (variance)
  df <- n - r - 1
  s2 <- SSerror / df
  # Valor esperado (cuando las observaciones son reales) 
  Xo <- c(1, Xo)  # añadir intercepto al vector de observaciones
  Yo <- sum(Xo * Beta)
  if (type == 1) {  # expected value
    LB <- Yo - sqrt(s2 * sum(Xo %*% W %*% Xo)) * qt((1 - alpha / 2), df)
    UB <- Yo + sqrt(s2 * sum(Xo %*% W %*% Xo)) * qt((1 - alpha / 2), df)
    Rows <- 'Valor esperado'
  }
  if (type == 2) {  # forecast
    LB <- Yo - sqrt(s2 * (1 + sum(Xo %*% W %*% Xo))) * qt((1 - alpha / 2), df)
    UB <- Yo + sqrt(s2 * (1 + sum(Xo %*% W %*% Xo))) * qt((1 - alpha / 2), df)
    Rows <- 'Forecast'
  }
  VarNames <- c('Respuesta', 'Lower Bound', 'Upper Bound')
  EVT <- data.frame(Yo = Yo, LB = LB, UB = UB, row.names = Rows)
  colnames(EVT) <- VarNames
  return(EVT)
}
