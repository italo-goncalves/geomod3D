#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// SPGP simulation
// [[Rcpp::export(.sparse_sim)]]
NumericVector sparse_sim(IntegerVector path, NumericVector nugget,
                         NumericVector w_, NumericMatrix Bi_,
                         NumericMatrix KMi_, NumericVector maxvar,
                         NumericMatrix K_, NumericVector d_,
                         NumericVector yTR, NumericVector vTR,
                         bool discount_noise, NumericVector Q_,
                         bool smooth, NumericVector randnum,
                         double reg = 1e-6){

  arma::mat K = arma::mat(K_.begin(), K_.nrow(), K_.ncol(), false);
  arma::mat Bi = arma::mat(Bi_.begin(), Bi_.nrow(), Bi_.ncol(), false);
  arma::mat KMi = arma::mat(KMi_.begin(), KMi_.nrow(), KMi_.ncol(), false);
  arma::colvec w = arma::colvec(w_.begin(), w_.length(), false);
  arma::colvec Q = arma::colvec(Q_.begin(), Q_.length(), false);
  arma::colvec yTR2 = arma::colvec(yTR.begin(), yTR.length(), false);
  arma::colvec nug = arma::colvec(nugget.begin(), nugget.length(), false);
  arma::colvec mvar = arma::colvec(maxvar.begin(), maxvar.length(), false);

  arma::mat I = arma::eye(size(Bi));
  arma::vec ysim = arma::vec(path.length());
  NumericVector vplus = d_ + vTR;

  path = path - 1;

  // sequential simulation
  for(int i = 0; i < path.length(); i++){
    arma::colvec k = K.row(path(i)).t();

    // prediction
    double mu = arma::as_scalar(sum(k % (Bi * w))); //+ yTR(path(i));
    // double v = maxvar - arma::as_scalar(sum(k % ((KMi - Bi) * k))) + vplus(path(i));
    double v = Q(path(i)) - arma::as_scalar(sum(k % ((KMi - Bi) * k))) + vplus(path(i));
    if (v < 0) v = 0; // to avoid rounding errors

    // simulation
    // ysim(path(i)) = rnorm(1, mu, sqrt(v))(0);
    ysim(path(i)) = randnum(path(i)) * sqrt(v) + mu;

    // update
    k = k / sqrt(d_(path(i)) + nug(path(i)));
    Bi = Bi * (I - (k * (k.t() * Bi)) /
      arma::as_scalar(1 + k.t() * (Bi * k)));
    k = k / sqrt(d_(path(i)) + nug(path(i)));
    w = w + ysim(path(i)) * k;
  }

  // smoothing (correction for varying quality of the sparse approximation)
  if(smooth){
    arma::mat dnew = 1 / (mvar - Q + reg);
    arma::mat Bnew = arma::inv(KMi) + K.t() * (arma::repmat(dnew, 1, K.n_cols) % K);
    Bnew = Bnew + arma::diagmat(arma::ones(K.n_cols) * reg); // regularization
    ysim = K * solve(Bnew, K.t() * (ysim % dnew));
  }

  // adding noise
  if(!discount_noise){
    arma::vec rand;
    rand.randn(ysim.n_elem);
    rand = rand % sqrt(nug);
    ysim = ysim + rand;
  }

  // end
  ysim = ysim + yTR2;
  return(wrap(ysim));
}

// [[Rcpp::export(.SPGP_CV)]]
List SPGP_CV(NumericVector w_, NumericMatrix Bi_, NumericVector y_,
                      NumericMatrix K_, NumericVector d_,
                      NumericVector yTR, NumericVector vTR){

  arma::mat K = arma::mat(K_.begin(), K_.nrow(), K_.ncol(), false);
  arma::mat Bi = arma::mat(Bi_.begin(), Bi_.nrow(), Bi_.ncol(), false);
  arma::colvec w = arma::colvec(w_.begin(), w_.length(), false);
  arma::colvec y = arma::colvec(y_.begin(), y_.length(), false);

  arma::mat I = arma::eye(size(Bi));
  arma::vec ymu = arma::vec(y_.length());
  arma::vec yvar = arma::vec(y_.length());
  NumericVector vplus = d_ + vTR;

  for(int i = 0; i < y_.length(); i++){
    arma::colvec k = K.row(i).t();

    // downdate
    arma::colvec k_tmp = k / sqrt(d_(i));
    arma::mat Bi_tmp = Bi * (I + (k_tmp * (k_tmp.t() * Bi)) /
      arma::as_scalar(1 - k_tmp.t() * (Bi * k_tmp)));
    k_tmp = k_tmp / sqrt(d_(i));
    arma::colvec w_tmp = w - y(i) * k_tmp;

    // prediction
    ymu(i) = arma::as_scalar(sum(k % (Bi_tmp * w_tmp))) + yTR(i);
    yvar(i) = arma::as_scalar(sum(k % (Bi_tmp * k))) + vplus(i);

  }

  return(List::create(Named("mean") = ymu,
                      Named("var") = yvar));
}
