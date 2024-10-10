mean_cov_corr_matrices <- function(X) {
#   %This function computes the Mean(M), the Variance - Covariance Matrix (VC)
#   %and the Correlation Matrix (R).
#   % The mean is computed by defining an n by 1 vector (1) that has equal angles
#   % with all the axes (observations)and the vector 1/(sqrt(n)) has unit
#   % lenght in the equal angle direction. Consider the vector yi = [x1i x2i
#   % ...xni](all values for one variable - column in X). The projection of yi
#   % on the vector 1 (by definition the proyection of x on y is (x'y/(y'y))*y)
#   % which es equal to (xi1 + x2i + ...xni)/n * 1 -> the mean vector. The mean
#   % is the proyection of the obeservation vector on the unit lenght vector of
#   % equal angles and the deviation from mena is the perperdicular distance of
#   % this two vectors. 
#   %The Variance - Covariance Matrix contains the variance 
#   %Skk = (1/(n-1)*sum((xji-mean(yi)^2)
#   
#   %in the diagonal and Sik = (1/(n-1)*sum((xji-mean(yi)*(xjk-mean(yk))
#   %outside the diagonal. The matrix is symetric.
#   %The sample correlation coefficient (linear association between two
#                                        %variables). Rik = Sik/(sqrt(Sii)sqrt(Skk)). 
#   % Is the standarization of the sample covariance provided by the square root of 
#   #sample variance of each variable.
#   %
#   %INPUT X a n by p matrix. n observation on p variables
#   %
#   %OUTPUT M Mean p by 1 vector
#   %       VC Variance - Covariance Matrix. p by p matrix 
#   %       R Correlation matrix. p by p matrix 
# 
# %--------------------------------------------------------------------------
  
  
  # Calcula la media M
  n <- nrow(X)
  p <- ncol(X)
  O <- matrix(1, n, 1)
  M <- (1/n) * t(X) %*% O
  
  # Matriz de desviaciones D
  Mr <- (1/n) * (O %*% t(O)) %*% as.matrix(X)
  D <- as.matrix(X) - Mr
  
  # Matriz de varianza-covarianza VC
  VC <- (1/(n-1)) * t(D) %*% D
  
  # Matriz de correlación R
  SD <- sqrt(diag(VC))
  SDinv <- solve(diag(SD))
  R <- SDinv %*% VC %*% SDinv
  
  return(list(M = M, VC = VC, R = R))
}