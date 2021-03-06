% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmnetALO.R
\name{glmnetALO}
\alias{glmnetALO}
\title{Compute ALO risk for glmnet object}
\usage{
glmnetALO(X, y, glm_obj, alpha, standardize, type.measure = c("mse",
  "mae", "deviance", "class"))
}
\arguments{
\item{X}{Input data matrix. Format should agree with that of the \code{glmnet} object.}

\item{y}{Response variable. Format should agree with that of the \code{glmnet} object.}

\item{glm_obj}{A fitted \code{glmnet} object.}

\item{alpha}{The elastic-net mixing parameter. Must agree with the value used to fit the \code{glmnet} object provided.}

\item{standardize}{Whether the \code{glmnet} object is fitted with x variable standardization.}

\item{type.measure}{Loss to use for CV risk estimation. Default to \code{"deviance"}.
\itemize{
  \item For \code{"gaussian"} family, one of \code{"mse"}, \code{"mae"}, and \code{"deviance"};
  \item For \code{"binomial"} family, one of \code{"mse"}, \code{"mae"}, \code{"class"}, and \code{"deviance"};
  \item For \code{"poisson"} family, one of \code{"mse"}, \code{"mae"}, and \code{"deviance"};
  \item For \code{"multinomial"} family, one of \code{"mse"}, \code{"mae"}, \code{"class"}, and \code{"deviance"}.
}}
}
\value{
An object with S3 class "glmnetALO".
\item{yALO}{The approximate leave-i-out fitted value.}
\item{beta}{The intercept and coefficients from the glmnet object. 
If \code{standardize = T} then unstandardized coefficients are returned. 
For multinomial models, coefficients returned are stacked column-wise, 
with the first K (number of classes) rows being the respective intercepts.}
\item{alom}{The mean ALO error.}
\item{alosd}{Estimate of standard error of \code{alom}.}
\item{aloup}{Upper curve = \code{alom + alosd}.}
\item{alolo}{Lower curve = \code{alom - alosd}.}
\item{nzeros}{Number of non-zero coefficients at each \code{lambda}.}
\item{alpha}{the values of alpha used in the fits.}
\item{lambda}{The values of \code{lambda} used in the fits.}
\item{lambda.min}{Value of \code{lambda} that gives minimum \code{alom}.}
\item{lambda.1se}{Largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
\item{family}{A text string indicating family of regression.}
\item{type.measure}{A text string indicating type of measure.}
}
\description{
This function computes the approximate leave-one-out risk estimation for the glmnet object provided.
}
\examples{
n = 300
p = 100
k = 60
beta = rnorm(p, 0, 1)
beta[-(1:k)] = 0
X = matrix(rnorm(n * p, 0, sqrt(1 / k)), n, p)
y = X \%*\% beta + rnorm(n, 0, 0.5)
fit = glmnet(X, y, family = "gaussian", alpha = 0.5, standardize = TRUE)
ALO_fit = glmnetALO(X, y, glm_obj = fit, alpha = 0.5, standardize = TRUE, type.measure = "mse")
}
\keyword{ALO}
