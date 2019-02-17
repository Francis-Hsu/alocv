#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat elnetALO(const arma::mat &XE, const arma::vec &yHat, const arma::mat &y, 
                   const double &lambda, const double &alpha, const bool &intercept) {
  arma::uword n = XE.n_rows; // # of samples
  arma::uword p = XE.n_cols; // # of active variables
  double sy = arma::as_scalar(arma::stddev(y, 1));
  double lScale = sy * n;
  
  arma::vec H;
  if (alpha == 1) {
    arma::mat R = chol(XE.t() * XE / lScale);
    arma::mat XER = solve(trimatl(R.t()), XE.t());
    H = sum(square(XER.t()), 1);
  } else {
    arma::mat hessR = (1 - alpha) * (lambda / std::pow(sy, 2)) * arma::eye<arma::mat>(p, p);
    if (intercept) {
      hessR(0, 0) = 0;
    }
    arma::mat F = inv_sympd(XE.t() * XE / lScale + hessR);
    H = sum((XE * F) % XE, 1);
  }
  arma::mat yALO = yHat + H % (yHat - y) / (lScale - H);
  
  
  return yALO;
}

arma::mat lognetALO(const arma::mat &XE, const arma::vec &yHat, const arma::mat &y, 
                    const double &lambda, const double &alpha, const bool &intercept) {
  arma::uword n = XE.n_rows; // # of samples
  arma::uword p = XE.n_cols; // # of active variables
  
  arma::vec eXb = exp(yHat);
  arma::mat hessR = (1 - alpha) * lambda * arma::eye<arma::mat>(p, p);
  if (intercept) {
    hessR(0, 0) = 0;
  }
  arma::mat F = inv_sympd(XE.t() * diagmat(eXb / square(1.0 + eXb)) * XE / n + hessR);
  arma::vec H = sum((XE * F) % XE, 1);
  arma::mat yALO = yHat + H % (eXb / (1.0 + eXb) - y) / (n - H % (eXb / square(1.0 + eXb)));
  
  return yALO;
}

arma::mat fishnetALO(const arma::mat &XE, const arma::vec &yHat, const arma::mat &y, 
                     const double &lambda, const double &alpha, const bool &intercept) {
  arma::uword n = XE.n_rows; // # of samples
  arma::uword p = XE.n_cols; // # of active variables
  
  arma::vec eXb = exp(yHat);
  arma::mat hessR = (1 - alpha) * lambda * arma::eye<arma::mat>(p, p);
  if (intercept) {
    hessR(0, 0) = 0;
  }
  arma::mat F = inv_sympd(XE.t() * diagmat(eXb) * XE / n + hessR);
  arma::vec H = sum((XE * F) % XE, 1);
  arma::mat yALO = yHat + H % (eXb - y) / (n - H % eXb);
  
  return yALO;
}

arma::mat multinetALO(const arma::mat &XESp, const arma::vec &yHat, const arma::mat &y, 
                      const double &lambda, const double &alpha, const bool &intercept) {
  // find out the dimension of X
  arma::uword K = y.n_cols; // # of classes
  arma::uword n = y.n_rows;
  arma::uword p = XESp.n_cols; // # of active variables, including intercepts (if any)
  
  // compute vector A(beta) and matrix D(beta)
  arma::vec A(n * K, arma::fill::zeros);
  arma::mat D(n * K, n * K, arma::fill::zeros);
  arma::vec tempA;
  for (arma::uword i = 0; i < n; ++i) {
    tempA = normalise(exp(yHat.rows(arma::span(i * K, i * K + K - 1))), 1); // normalise to have unit 1-norm
    A(arma::span(i * K, i * K + K - 1)) = tempA;
    D(arma::span(i * K, i * K + K - 1), arma::span(i * K, i * K + K - 1)) = diagmat(tempA) - tempA * tempA.t();
  }
  
  // compute hessR
  arma::mat hessR = n * lambda * (1 - alpha) * arma::eye<arma::mat>(p, p);
  if (intercept) {
    hessR(arma::span(0, K - 1), arma::span(0, K - 1)) = arma::zeros<arma::mat>(K, K); // no penalty for intercepts
  }
  
  // compute matrix K(beta) and its inverse
  arma::mat K_inv = pinv(XESp.t() * D * XESp + hessR);
  
  // do leave-i-out prediction
  arma::mat aloUpdate(K, n);
  arma::mat X_i, D_i, XKX, XKXAmy;
  for (arma::uword i = 0; i < n; i++) {
    // find the X_i and y_i
    X_i = XESp.rows(arma::span(i * K, i * K + K - 1));
    
    // find D_i
    D_i = D(arma::span(i * K, i * K + K - 1), arma::span(i * K, i * K + K - 1));
    
    // compute XKX
    XKX = X_i * K_inv * X_i.t();
    XKXAmy = XKX * (A(arma::span(i * K, i * K + K - 1)) - arma::trans(y.row(i)));
    
    // compute the leave-i-out update
    aloUpdate.col(i) = XKXAmy - XKX * pinv(XKX - pinv(D_i)) * XKXAmy;
  }
  arma::mat yALO = yHat + vectorise(aloUpdate);
  
  return(yALO);
}

//[[Rcpp::export]]
arma::mat glmnetALODirect(const arma::mat &X, const arma::mat &y, 
                          const arma::sp_mat &beta, const arma::vec &lambda, const double &alpha,
                          const arma::field<arma::uvec> &addList, 
                          const arma::field<arma::uvec> &dropList,
                          const int family, const bool intercept) {
  // function pointer for polymorphism
  arma::mat (*netPtr)(const arma::mat &, const arma::vec &, const arma::mat &, const double &, const double &, const bool &) = NULL;
  if (family == 0) {
    netPtr = &elnetALO;
  } else if (family == 1) {
    netPtr = &lognetALO;
  } else if (family == 2) {
    netPtr = &fishnetALO;
  } else if (family == 3) {
    netPtr = &multinetALO;
  }
  
  // preamble
  arma::uvec currActiveID = addList(0);
  arma::mat yALO(X.n_rows, beta.n_cols, arma::fill::zeros); // matrix storing ALO estimates
  arma::vec yHat = X * beta.col(0);
  arma::mat XE = X.cols(currActiveID);
  yALO.col(0) = netPtr(XE, yHat, y, lambda(0), alpha, intercept);
  
  // update through all lambdas
  arma::uvec toAdd;
  arma::uvec toDrop;
  arma::uvec keepID;
  arma::uvec dropID;
  for (arma::uword i = 1; i < beta.n_cols; i++) {
    // checkUserInterrupt();
    toAdd = addList(i); // variables to add
    toDrop = dropList(i - 1); // variables to drop
    
    if (currActiveID.n_elem > 0) {
      // remove variables
      if (toDrop.n_elem > 0) {
        keepID = arma::regspace<arma::uvec>(0, currActiveID.n_elem - 1);
        dropID.zeros(toDrop.n_elem);
        for (arma::uword j = 0; j < toDrop.n_elem; j++) {
          dropID(j) = arma::as_scalar(arma::find(currActiveID == toDrop(j), 1, "first"));
        }
        dropID = sort(dropID, "descend");
        for (arma::uword j = 0; j < toDrop.n_elem; j++) {
          currActiveID.shed_row(dropID(j));
          keepID.shed_row(dropID(j));
        }
        XE = XE.cols(keepID);
      }
      // add variables
      if (toAdd.n_elem > 0) {
        currActiveID = join_cols(currActiveID, toAdd);
        XE = join_rows(XE, X.cols(toAdd));
      }
    } else {
      currActiveID = toAdd;
      XE = X.cols(toAdd);
    }
    yHat = X * beta.col(i);
    yALO.col(i) = netPtr(XE, yHat, y, lambda(i), alpha, intercept);
  }
  
  return yALO;
}

// [[Rcpp::export]]
arma::mat multinetExpand(const arma::mat &X, const arma::uword K) {
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  arma::mat XSp(n * K, p * K, arma::fill::zeros);
  arma::uvec rowID, colID;
  for (arma::uword i = 0; i < K; i++) {
    rowID = arma::regspace<arma::uvec>(0, K, n * K - K) + i;
    colID = arma::regspace<arma::uvec>(p * i, p * i + p - 1);
    XSp.submat(rowID, colID) = X;
  }
  
  return XSp;
}