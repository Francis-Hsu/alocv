#' Plot the ALO curve produced by \code{glmnetALO}
#'
#' Plots the ALO curve, and upper and lower standard deviation curves, as a function of the \code{lambda} values used.
#' @param x Fitted \code{\link{glmnetALO}} object.
#' @param ... Other graphical parameters to plot.
#' @importFrom graphics abline axis points segments
#' @export
plot.glmnetALO = function(x, ...) {
  # extract useful stuffs
  obj_measure = x$type.measure
  obj_family = x$family
  obj_riskm = x$alom
  obj_riskup = x$aloup
  obj_risklo = x$alolo
  obj_lambda = x$lambda
  obj_lambda.min = x$lambda.min
  obj_lambda.1se = x$lambda.1se
  obj_nz = x$nzeros
  
  log_param = log(obj_lambda)
  
  xlab = "log(Lambda)"
  ylab = switch(obj_measure,
                mse = "Mean-Squared Error",
                mae = "Mean Absolute Error",
                deviance = "Deviance",
                class = "Class")
  plot.args = list(x = log(obj_lambda),
                   y = obj_riskm,
                   ylim = range(obj_riskup, obj_risklo), 
                   xlab = xlab, 
                   ylab = ylab, 
                   type = "n")
  new.args = list(...)
  if (length(new.args)){ 
    plot.args[names(new.args)] = new.args
  }
  do.call("plot", plot.args)
  
  # plot error bars
  barw = diff(range(log_param)) * 0.01
  segments(log_param, obj_riskup, log_param, obj_risklo, col = "darkgrey")
  segments(log_param - barw, obj_riskup, log_param + barw, obj_riskup, col = "darkgrey")
  segments(log_param - barw, obj_risklo, log_param + barw, obj_risklo, col = "darkgrey")
  
  # plot risks
  points(log_param, obj_riskm, pch = 20, col = "blue")
  axis(side = 3, at = log_param, labels = paste(obj_nz), tick = FALSE, line = 0)
  abline(v = log(obj_lambda.min), lty = 3)
  abline(v = log(obj_lambda.1se), lty = 3)
  invisible()
}