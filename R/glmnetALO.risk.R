glmnetALO.risk = function(yalo, y, family, type.measure) {
  if (family == 0) {
    errmat = sweep(yalo, 1, y)
    aloraw = switch(type.measure,
                    mse = errmat^2,
                    deviance = errmat^2,
                    mae = abs(errmat))
  } else if (family == 1) {
    y = as.factor(y)
    nc = as.integer(length(table(y)))
    y = diag(nc)[as.numeric(y), ]
    predmat = 1 / (1 + exp(-yalo))
    prob_min = 1e-05
    prob_max = 1 - prob_min
    aloraw = switch(type.measure,
                    mse = (y[, 1] - (1 - predmat))^2 + (y[, 2] - predmat)^2,
                    mae = abs(y[, 1] - (1 - predmat)) + abs(y[, 2] - predmat),
                    class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <=0.5),
                    deviance = {
                      predmat = pmin(pmax(predmat, prob_min), prob_max)
                      lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                      ly = log(y)
                      ly[y == 0] = 0
                      ly = drop((y * ly) %*% c(1, 1))
                      2 * (ly - lp)
                    })
  } else if (family == 2) {
    errmat = sweep(exp(yalo), 1, y)
    aloraw = switch(type.measure,
                    mse = errmat^2,
                    mae = abs(errmat))
  } else if (family == 3) {
    m = ncol(yalo)
    n = nrow(y)
    K = ncol(y)
    errmat = matrix(0, n * K, m)
    for (i in 1:m) {
      predmat = exp(matrix(yalo[, i], nrow(y), ncol(y), byrow = T))
      predmat = predmat / rowSums(predmat)
      errmat[, i] = c(y) - c(predmat)
    }
    aloraw = switch(type.measure,
                    mse = K * errmat^2,
                    mae = K * abs(errmat))
  }
  alom = colMeans(aloraw)
  alosd = sqrt(colMeans(scale(aloraw, alom, FALSE)^2, na.rm = TRUE) / (length(y) - 1))
  
  return(list(alom = alom, alosd = alosd))
}