#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]

// lambda - penalty
// X, y - data and response
// w - svm coefficent
// b - svm intercept
// b_scaling - intercept scaling factor
// tol - tolerance for detecting singularities
//[[Rcpp::export]]
vec svmALO_LIBLINEAR(const double &lambda, const mat &X, const vec &y, const vec &w, const double &b,
                     const double &b_scaling, const double &tol) {
  // augment the data and weight matrices with bias
  mat B = b_scaling * ones<vec>(X.n_rows);
  mat XAug = join_rows(X, B);
  vec wAug(w.n_elem + 1);
  wAug(span(0, w.n_elem - 1)) = w;
  wAug(w.n_elem) = b;

  vec yHat = XAug * wAug;
  vec yyHat = y % yHat;

  // identify singularities
  uvec V = find(abs(1 - yyHat) < tol);
  uvec S = find(abs(1 - yyHat) >= tol);

  mat I_p = eye<mat>(XAug.n_cols, XAug.n_cols);
  mat XV = XAug.rows(V);
  mat inv_XXV = inv_sympd(XV * XV.t());

  uvec gIdx = intersect(find(yyHat < 1.0), S);
  vec y_gIdx = y.elem(gIdx);
  mat X_gIdx = XAug.rows(gIdx);

  // containers for a and g
  vec a = zeros<vec>(XAug.n_rows);
  vec g = zeros<vec>(XAug.n_rows);

  // compute a and g for S
  mat Xa_s = XAug.rows(S) * (I_p - XV.t() * inv_XXV * XV) * trans(XAug.rows(S)) / lambda;
  a.elem(S) = diagvec(Xa_s);
  g.elem(gIdx) = -y.elem(gIdx);

  // compute a and g for V
  vec l_beta = lambda * wAug;
  vec sum_yX = trans(sum(diagmat(y_gIdx) * X_gIdx, 0));
  a.elem(V) = 1 / (lambda * diagvec(inv_XXV));
  g.elem(V) = inv_XXV * XV * (l_beta - sum_yX);

  return yHat + a % g;
}

//[[Rcpp::export]]
vec svmALO(const mat &X, const vec &y, const vec &w, const double &b, const double &lambda,
           const double &tol) {
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
vec svmKerALO(const mat &K, const vec &y, const vec &alpha, const double &rho, const double &lambda,
              const double &tol) {
  arma::uword N = y.n_elem;

  // augment the data and weight matrices with offset
  arma::mat Kinv = arma::inv(K); // cannot use inv_sympd because of sigmoid kernel

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
vec nuSvcKerALO(const mat &K, const vec &y, const vec &alpha, const double &rho, const double &nu, const double &r,
                const double &tol) {
  arma::uword N = y.n_elem;
  
  // augment the data and weight matrices with offset
  arma::mat Kinv = arma::inv(K); // cannot use inv_sympd because of sigmoid kernel
  
  arma::vec yHat = K * alpha - rho;
  arma::vec yyHat = y % yHat;
  
  // identify singularities
  arma::uvec V = arma::find(arma::abs(r - yyHat) < tol);
  arma::uvec S = arma::find(arma::abs(r - yyHat) >= tol);
  
  // useful matrices
  arma::mat I_n = arma::eye<mat>(N, N);
  arma::mat KV = K.cols(V);
  arma::mat KS = K.cols(S);
  arma::mat K1 = arma::inv(KV.t() * Kinv * KV);
  
  arma::uvec gID = arma::intersect(arma::find(yyHat < r), S);
  arma::mat yK_g = K.cols(gID);
  yK_g.each_row() %= arma::trans(y.elem(gID));
  
  // containers for a and g
  arma::vec a = zeros<vec>(N);
  arma::vec g = zeros<vec>(N);
  
  // compute a and g for S
  arma::mat Ka_s = KS.t() * Kinv * (I_n - KV * K1 * KV.t() * Kinv) * KS * (nu * r);
  a.elem(S) = arma::diagvec(Ka_s);
  g.elem(gID) = -y.elem(gID) / N;
  a.print();
  // compute a and g for V
  arma::vec gradR = K * alpha;
  arma::vec sum_yK = arma::sum(yK_g, 1) / N;
  a.elem(V) = (nu * r) / arma::diagvec(K1);
  g.elem(V) = arma::inv(KV.t() * KV) * KV.t() * (sum_yK  - gradR) / (nu * r);
  
  return yHat + a % g;
}


//[[Rcpp::export]]
vec svrKerALO(const mat &K, const vec &y, const vec &alpha, const double &rho, const double &lambda, const double &epsilon,
              const double &tol) {
  arma::uword N = y.n_elem;
  
  // augment the data and weight matrices with offset
  arma::mat Kinv = arma::inv(K); // cannot use ALO_risk because of sigmoid kernel
  
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
  g.elem(V) = arma::inv(KV.t() * KV) * KV.t() * (gradR - sum_sgnK);
  
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
arma::mat sigmoidKer(const arma::mat &X, const double &gamma, const double &coef0) {
  // compute the kernel matrix
  arma::mat K = arma::tanh(gamma * X * X.t() + coef0);

  return K;
}

//[[Rcpp::export]]
arma::mat polynomialKer(const arma::mat &X, const double &gamma, const double &coef0, const int &degree) {
  // compute the kernel matrix
  arma::mat K = arma::pow(gamma * X * X.t() + coef0, degree);

  return K;
}

//[[Rcpp::export]]
arma::mat gaussianKerApprox(const arma::mat &X, const double &gamma) {
  arma::mat K = gaussianKer(X, gamma);

  // SVD to find matrix square root
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, K);
  s = arma::clamp(s, 1e-12, s.max());
  U.each_row() /= arma::sqrt(s).t();

  return K * U * V;
}

//[[Rcpp::export]]
arma::mat sigmoidKerApprox(const arma::mat &X, const double &gamma, const double &coef0) {
  // compute the kernel matrix
  arma::mat K = sigmoidKer(X, gamma, coef0);

  // SVD to find matrix square root
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, K);
  s = arma::clamp(s, 1e-12, s.max());
  U.each_row() /= arma::sqrt(s).t();

  return K * U * V;
}

//[[Rcpp::export]]
arma::mat polynomialKerApprox(const arma::mat &X, const double &gamma, const double &coef0, const int &degree) {
  // compute the kernel matrix
  arma::mat K = polynomialKer(X, gamma, coef0, degree);

  // SVD to find matrix square root
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, K);
  s = arma::clamp(s, 1e-12, s.max());
  U.each_row() /= arma::sqrt(s).t();

  return K * U * V;
}
