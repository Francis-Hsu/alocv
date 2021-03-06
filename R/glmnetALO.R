#' Compute ALO risk for glmnet object
#'
#' This function computes the approximate leave-one-out risk estimation for the glmnet object provided.
#' @param X Input data matrix. Format should agree with that of the \code{glmnet} object.
#' @param y Response variable. Format should agree with that of the \code{glmnet} object.
#' @param glm_obj A fitted \code{glmnet} object.
#' @param alpha The elastic-net mixing parameter. Must agree with the value used to fit the \code{glmnet} object provided.
#' @param standardize Whether the \code{glmnet} object is fitted with x variable standardization.
#' @param type.measure Loss to use for CV risk estimation. Default to \code{"deviance"}.
#' \itemize{
#'   \item For \code{"gaussian"} family, one of \code{"mse"}, \code{"mae"}, and \code{"deviance"};
#'   \item For \code{"binomial"} family, one of \code{"mse"}, \code{"mae"}, \code{"class"}, and \code{"deviance"};
#'   \item For \code{"poisson"} family, one of \code{"mse"}, \code{"mae"}, and \code{"deviance"};
#'   \item For \code{"multinomial"} family, one of \code{"mse"}, \code{"mae"}, \code{"class"}, and \code{"deviance"}.
#' }
#' @return An object with S3 class "glmnetALO".
#' \item{yALO}{The approximate leave-i-out fitted value.}
#' \item{beta}{The intercept and coefficients from the glmnet object. 
#' If \code{standardize = T} then unstandardized coefficients are returned. 
#' For multinomial models, coefficients returned are stacked column-wise, 
#' with the first K (number of classes) rows being the respective intercepts.}
#' \item{alom}{The mean ALO error.}
#' \item{alosd}{Estimate of standard error of \code{alom}.}
#' \item{aloup}{Upper curve = \code{alom + alosd}.}
#' \item{alolo}{Lower curve = \code{alom - alosd}.}
#' \item{nzeros}{Number of non-zero coefficients at each \code{lambda}.}
#' \item{alpha}{the values of alpha used in the fits.}
#' \item{lambda}{The values of \code{lambda} used in the fits.}
#' \item{lambda.min}{Value of \code{lambda} that gives minimum \code{alom}.}
#' \item{lambda.1se}{Largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
#' \item{family}{A text string indicating family of regression.}
#' \item{type.measure}{A text string indicating type of measure.}
#' @keywords ALO
#' @export
#' @examples
#' n = 300
#' p = 100
#' k = 60
#' beta = rnorm(p, 0, 1)
#' beta[-(1:k)] = 0
#' X = matrix(rnorm(n * p, 0, sqrt(1 / k)), n, p)
#' y = X %*% beta + rnorm(n, 0, 0.5)
#' fit = glmnet(X, y, family = "gaussian", alpha = 0.5, standardize = TRUE)
#' ALO_fit = glmnetALO(X, y, glm_obj = fit, alpha = 0.5, standardize = TRUE, type.measure = "mse")

glmnetALO = function(X, y, glm_obj, alpha, standardize, type.measure = c("mse", "mae", "deviance", "class")) {
  # extract useful stuffs from the glmnet object
  glm_family = class(glm_obj)[1]
  glm_lambda = glm_obj$lambda # assume at least 2 lambdas
  glm_beta = glm_obj$beta # coefficients
  glm_a0 = glm_obj$a0 # intercept(s)
  nz = glm_obj$df # number of nonzero coefs
  intercept_flag = any(glm_a0 != 0) # check if the model has intercept
  
  # response type of the glmnet object
  # support "gaussian", "binomial", "poisson" and "multinomial" for now
  family_id = switch(glm_family,
                      elnet = 0,
                      lognet = 1,
                      fishnet = 2,
                      multnet = 3,
                      stop("Unsupported regression family, should be one of
                           \"gaussian\", \"binomial\", \"poisson\", or \"multinomial\""))
  
  # validate type of measure to use
  if (missing(type.measure)) {
    type.measure = "deviance"
  } else {
    type.measure = match.arg(type.measure)
  }
  
  # Preprocessing
  # coerce input format
  X = as.matrix(X)
  if(family_id == 1 || family_id == 3) {
    if(!is.matrix(y)) {
      y = as.factor(y)
      nc = as.integer(length(table(y)))
      y = diag(nc)[as.numeric(y), ]
    }
  } else {
    y = as.matrix(y)
  }
  if(family_id == 1) {
    y = y[, 2, drop = F] # glmnet treated the second column as the target class
  }
  n = nrow(X) # number of observations
  K = ncol(y) # number of classes, for multinomial regression
  
  # standardize input and restandardize coefficients
  if (standardize) {
    glm_rescale = glmnetALO.rescale(X, glm_a0, glm_beta, intercept_flag, family_id)
    X = glm_rescale$Xs
    glm_a0 = glm_rescale$a0s
    glm_beta = glm_rescale$betas
  }
  
  # prepare covariates and coefficients
  if (family_id == 3) {
    # convert factor to indicator matrix
    X = multinetExpand(X, K)
    glm_beta = do.call(rbind, glm_beta)
    if (intercept_flag) {
      X = cbind(do.call(rbind, replicate(nrow(y), diag(K), simplify = FALSE)), X)
      glm_beta = rbind(glm_a0, glm_beta)
    }
  } else {
    if (intercept_flag) {
      X = cbind(1, X)
      glm_beta = rbind(glm_a0, glm_beta)
    }
  }
  # glm_beta = unname(glm_beta) # drop the names
  
  # identify active set, find variables to add/remove
  M = length(glm_lambda)
  glm_active_set = lapply(1:M, function(i) { which(abs(glm_beta[, i]) >= 1e-8) - 1 })
  glm_add_idx = lapply(2:M, function(i) {
    setdiff(glm_active_set[[i]], glm_active_set[[i - 1]])
  })
  glm_add_idx = append(glm_add_idx, list(glm_active_set[[1]]), 0)
  glm_rmv_idx = lapply(2:M, function(i) {
    setdiff(glm_active_set[[i - 1]], glm_active_set[[i]])
  })
  
  # compute ALO and corresponding risk
  y_alo = glmnetALODirect(X, y, glm_beta, glm_lambda, alpha, 
                          glm_add_idx, glm_rmv_idx, 
                          family_id, intercept_flag)
  alo_risk = glmnetALO.risk(y_alo, y, family_id, type.measure)
  alom = alo_risk$alom
  alosd = alo_risk$alosd
  
  # find optimum values of lambda
  idmin = which.min(alo_risk$alom)
  lambda.min = glm_lambda[idmin]
  semin = alom[idmin] + alosd[idmin]
  id1se = which(alom <= semin)
  lambda.1se = max(glm_lambda[id1se], na.rm = TRUE)
  
  # create and return S3 glmnetALO object
  # do we need S4?
  rtn_obj = list(yALO = y_alo, beta = glm_beta,
                 alom = alom, alosd = alosd, aloup = alom + alosd, alolo = alom - alosd,
                 nzeros = nz, alpha = alpha,
                 lambda = glm_lambda, lambda.min = lambda.min, lambda.1se = lambda.1se, 
                 family = glm_family, type.measure = type.measure)
  attr(rtn_obj, "class") = "glmnetALO"
  return(rtn_obj)
}