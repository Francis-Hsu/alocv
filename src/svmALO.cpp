#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::vec svcALO(const arma::mat &X, const arma::vec &y, const arma::vec &w, 
                 const double &b, const double &lambda, const double &tol) {
  arma::uword N = X.n_rows;
  arma::uword P = X.n_cols;
  arma::vec yHat = X * w + b;
  arma::vec yyHat = y % yHat;

  // identify singularities
  arma::uvec V = arma::find(arma::abs(1 - yyHat) < tol);
  arma::uvec S = arma::find(arma::abs(1 - yyHat) >= tol);

  // useful matrices
  arma::mat I_p = arma::eye<mat>(P, P);
  arma::mat XV = X.rows(V);
  arma::mat XS = X.rows(S);
  arma::mat inv_XXV = arma::inv(XV * XV.t());

  arma::uvec gID = arma::intersect(arma::find(yyHat < 1.0), S);
  arma::mat yX_g = X.rows(gID);
  yX_g.each_col() %= y.elem(gID);

  // containers for a and g
  arma::vec a = zeros<vec>(N);
  arma::vec g = zeros<vec>(N);

  // compute a and g for S
  arma::mat Xa_s = XS * (I_p - XV.t() * inv_XXV * XV) * XS.t() / lambda;
  a.elem(S) = arma::diagvec(Xa_s);
  g.elem(gID) = -y.elem(gID);

  // compute a and g for V
  arma::vec gradR = lambda * w;
  arma::vec sum_yX = trans(arma::sum(yX_g, 0));
  a.elem(V) = 1 / (lambda * arma::diagvec(inv_XXV));
  g.elem(V) = inv_XXV * XV * (sum_yX - gradR);

  return yHat + a % g;
}

//[[Rcpp::export]]
arma::vec svcKerALO(const arma::mat &K, const arma::vec &y, const arma::vec &alpha, 
                    const double &rho, const double &lambda, const double &tol) {
  arma::uword N = y.n_elem;

  // augment the data and weight matrices with offset
  arma::mat Kinv = arma::inv_sympd(K);

  arma::vec yHat = K * alpha - rho;
  arma::vec yyHat = y % yHat;

  // identify singularities
  arma::uvec V = arma::find(arma::abs(1 - yyHat) < tol);
  arma::uvec S = arma::find(arma::abs(1 - yyHat) >= tol);

  // useful matrices
  arma::mat I_n = arma::eye<mat>(N, N);
  arma::mat KV = K.cols(V);
  arma::mat KS = K.cols(S);
  arma::mat K1 = arma::inv(KV.t() * Kinv * KV);

  arma::uvec gID = arma::intersect(arma::find(yyHat < 1.0), S);
  arma::mat yK_g = K.cols(gID);
  yK_g.each_row() %= arma::trans(y.elem(gID));

  // containers for a and g
  arma::vec a = zeros<vec>(N);
  arma::vec g = zeros<vec>(N);

  // compute a and g for S
  arma::mat Ka_s = KS.t() * Kinv * (I_n - KV * K1 * KV.t() * Kinv) * KS / lambda;
  a.elem(S) = arma::diagvec(Ka_s);
  g.elem(gID) = -y.elem(gID);

  // compute a and g for V
  arma::vec gradR = lambda * K * alpha;
  arma::vec sum_yK = arma::sum(yK_g, 1);
  a.elem(V) = 1 / (lambda * arma::diagvec(K1));
  g.elem(V) = arma::inv(KV.t() * KV) * KV.t() * (sum_yK - gradR);

  return yHat + a % g;
}

//[[Rcpp::export]]
arma::vec svrKerALO(const arma::mat &K, const arma::vec &y, const arma::vec &alpha, 
                    const double &rho, const double &lambda, const double &epsilon, const double &tol) {
  arma::uword N = y.n_elem;
  
  // augment the data and weight matrices with offset
  arma::mat Kinv = arma::inv_sympd(K); // not gonna work for sigmoid kernel
  
  arma::vec yHat = K * alpha - rho;
  arma::vec absResid = arma::abs(y - yHat);
  
  // identify singularities
  arma::uvec V = arma::find(arma::abs(absResid - epsilon) < tol);
  arma::uvec S = arma::find(arma::abs(absResid - epsilon) >= tol);
  
  // useful matrices
  arma::mat I_n = arma::eye<mat>(N, N);
  arma::mat KV = K.cols(V);
  arma::mat KS = K.cols(S);
  arma::mat K1 = arma::inv(KV.t() * Kinv * KV);
  
  arma::uvec gID = arma::intersect(arma::find(absResid >= epsilon), S);
  arma::mat sgnK_g = K.cols(gID);
  sgnK_g.each_row() %= arma::trans(arma::sign(yHat.elem(gID)));
  
  // containers for a and g
  arma::vec a = zeros<vec>(N);
  arma::vec g = zeros<vec>(N);
  
  // compute a and g for S
  arma::mat Ka_s = KS.t() * Kinv * (I_n - KV * K1 * KV.t() * Kinv) * KS / lambda;
  a.elem(S) = arma::diagvec(Ka_s);
  g.elem(gID) = -arma::sign(yHat.elem(gID));
  
  // compute a and g for V
  arma::vec gradR = lambda * K * alpha;
  arma::vec sum_sgnK = arma::sum(sgnK_g, 1);
  a.elem(V) = 1 / (lambda * arma::diagvec(K1));
  g.elem(V) = arma::inv(KV.t() * KV) * KV.t() * (sum_sgnK - gradR);
  
  return yHat + a % g;
}

//[[Rcpp::export]]
arma::mat gaussianKer(const arma::mat &X, const double &gamma) {
  // compute the Euclidean distance matrix
  arma::vec rowSumSq = arma::sum(arma::square(X), 1);
  arma::mat D = -2 * X * X.t();
  D.each_row() += rowSumSq.t();
  D.each_col() += rowSumSq;

  return arma::exp(-gamma * D);
}

//[[Rcpp::export]]
arma::mat polynomialKer(const arma::mat &X, const double &gamma, const double &coef0, const int &degree) {
  // compute the kernel matrix
  arma::mat K = arma::pow(gamma * X * X.t() + coef0, degree);

  return K;
}

// arma::mat sigmoidKer(const arma::mat &X, const double &gamma, const double &coef0) {
//   // compute the kernel matrix
//   arma::mat K = arma::tanh(gamma * X * X.t() + coef0);
// 
//   return K;
// }