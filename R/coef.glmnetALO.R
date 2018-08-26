#' Extract the coefficients from a \code{glmnetALO} object
#'
#' This function extracts coefficients from a \code{glmnetALO} object for the optimal value chosen for \code{lambda}.
#' @param obj Fitted \code{\link{glmnetALO}} object.
#' @param s Value of the penalty parameter \code{lambda} at which coefficients are required. 
#' Default is the value \code{s = "lambda.1se"}. Alternatively \code{s = "lambda.min"} can be used. 
#' If \code{s} is numeric, it is taken as the value of \code{lambda} to be used.
#' @return An matrix of coefficients, nothing will be returned if no matching \code{s} is found.
#' @export
coef.glmnetALO = function(obj, s = c("lambda.1se", "lambda.min")) {
  obj_lambda = obj$lambda
  
  if(is.numeric(s)) {
    which = match(s, obj_lambda, FALSE)
  } else if (is.character(s)) {
    s = match.arg(s)
    which = match(obj[[s]], obj_lambda, FALSE)
  } else {
    stop("Invalid form for s.")
  }
  if (!which) {
    stop("No matching s was found.")
  }
  
  return(obj$beta[, which, drop = F])
}