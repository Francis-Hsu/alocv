# estimate risks for fitted ALO value
# most are taken direcly from glmnet package (2.0-16), with minor modifications
glmnetALO.risk = function(yalo, y, family, type.measure) {
  if (family == 0) {# linear
    errmat = sweep(yalo, 1, y)
    aloraw = switch(type.measure,
                    mse = errmat^2,
                    deviance = errmat^2,
                    mae = abs(errmat))
  } else if (family == 1) {# logistic
    # convert to indicator matrix
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
  } else if (family == 2) {# Poisson
    aloraw = switch(type.measure,
                    mse = (sweep(exp(yalo), 1, y))^2,
                    mae = abs(sweep(exp(yalo), 1, y)),
                    deviance = {
                      deveta = c(y) * yalo - exp(yalo)
                      devy = y * log(y) - y
                      devy[y == 0] = 0
                      2 * sweep(-deveta, 1, devy, "+")
                    })
  } else if (family == 3) {# multinomial
    m = ncol(yalo)
    n = nrow(y)
    K = ncol(y)
    prob_min = 1e-05
    prob_max = 1 - prob_min
    predmat = array(0, dim = c(n, K, m))
    for (i in 1:m) {
      predmat[, , i] = exp(matrix(yalo[, i], n, K, byrow = T))
      predmat[, , i] = predmat[, , i] / rowSums(predmat[, , i])
    }
    bigY = array(y, dim(predmat))
    aloraw = switch(type.measure,
                    mse = apply((bigY - predmat)^2, c(1, 3), sum),
                    mae = apply(abs(bigY - predmat), c(1, 3), sum),
                    deviance = {
                      predmat = pmin(pmax(predmat, prob_min), prob_max)
                      lp = bigY * log(predmat)
                      ly = bigY * log(bigY)
                      ly[bigY == 0] = 0
                      apply(2 * (ly - lp), c(1, 3), sum)
                    }, 
                    class = {
                      classid = c(as.numeric(apply(predmat, c(1, 3), which.max)))
                      yperm = matrix(aperm(bigY, c(1, 3, 2)), ncol = K)
                      matrix(1 - yperm[cbind(seq(classid), classid)], ncol = m)
                    })
  }
  alom = colMeans(aloraw)
  alosd = sqrt(colMeans(scale(aloraw, alom, FALSE)^2, na.rm = TRUE) / (length(y) - 1))
  
  return(list(alom = alom, alosd = alosd))
}