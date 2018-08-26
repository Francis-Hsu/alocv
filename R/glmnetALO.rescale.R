glmnetALO.rescale = function(X, a0, beta, intercept, family) {
  mean_X = colMeans(X)
  sd_X = sqrt(colSums(X^2 / n) - colSums(X / n)^2)
  X = scale(X, center = intercept, scale = sd_X) # no centering if no intercept
  
  if (family == 3) {
    beta = lapply(beta, "*", sd_X)
    for (i in 1:nrow(a0)) {
      a0[i, ] = as.vector(a0[i, ] + mean_X %*% beta[[i]])
    }
  } else {
    beta = beta * sd_X
    a0 = as.vector(a0 + mean_X %*% beta)
  }
  
  return(list(Xs = X, a0s = a0, betas = beta))
}