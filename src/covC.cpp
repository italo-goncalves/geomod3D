#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export(.vectorized_pdist)]]
NumericMatrix vectorized_pdist(NumericMatrix Ar, NumericMatrix Br) {
	// code from http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
    int m = Ar.nrow(),
        n = Br.nrow(),
        k = Ar.ncol();
    arma::mat A = arma::mat(Ar.begin(), m, k, false);
    arma::mat B = arma::mat(Br.begin(), n, k, false);

    arma::colvec An =  sum(square(A),1);
    arma::colvec Bn =  sum(square(B),1);

    arma::mat C = -2 * (A * B.t());
    C.each_col() += An;
    C.each_row() += Bn.t();

    // avoiding NaNs - code above can give negative squared distances
    C.elem( arma::find(C < 0) ).zeros();

    return wrap(sqrt(C));
}

// [[Rcpp::export(.covd1_gaussian)]]
NumericMatrix covd1_gaussian(NumericMatrix u, NumericMatrix v,
                             NumericMatrix dir1, NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec dvec = D1.row(j).t();

                        double K = arma::as_scalar(exp(-3 * sum(square(B *
                                                   (uvec - vvec)))));
                        arma::colvec h = uvec - vvec;
                        arma::colvec grad = 6 * K * AtA * h;

                        C(i, j) = arma::as_scalar(grad.t() * dvec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd2_gaussian)]]
NumericMatrix covd2_gaussian(NumericMatrix u, NumericMatrix v,
                             NumericMatrix dir1, NumericMatrix dir2,
                             NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);
        arma::mat D2 = arma::mat(dir2.begin(), dir2.nrow(), dir2.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec d1vec = D1.row(i).t();
                        arma::colvec d2vec = D2.row(j).t();

                        double K = arma::as_scalar(exp(-3 * sum(square(B *
                                                   (uvec - vvec)))));
                        arma::colvec h = AtA * (uvec - vvec);
                        arma::mat hess = K * (6 * AtA - 36 * h * h.t());

                        C(i, j) = arma::as_scalar(d1vec.t() * hess * d2vec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd1_cubic)]]
NumericMatrix covd1_cubic(NumericMatrix u, NumericMatrix v,
                             NumericMatrix dir1, NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec dvec = D1.row(j).t();

                        double d = arma::as_scalar(sqrt(sum(square(B *
                                                   (uvec - vvec)))));
                        d += 1e-12;
                        if(d > 1){
                                C(i, j) = 0;
                        }
                        else{
                                arma::colvec h = uvec - vvec;
                                arma::colvec grad = (14 - 105/4 * d +
                                        35/2 * pow(d,3) - 21/4 * pow(d,5)) *
                                        AtA * h;

                                C(i, j) = arma::as_scalar(grad.t() * dvec);
                        }
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd2_cubic)]]
NumericMatrix covd2_cubic(NumericMatrix u, NumericMatrix v,
                             NumericMatrix dir1, NumericMatrix dir2,
                             NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);
        arma::mat D2 = arma::mat(dir2.begin(), dir2.nrow(), dir2.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec d1vec = D1.row(i).t();
                        arma::colvec d2vec = D2.row(j).t();

                        double d = arma::as_scalar(sqrt(sum(square(B *
                                                   (uvec - vvec)))));
                        d += 1e-12;
                        if(d > 1){
                                C(i, j) = 0;
                        }
                        else{
                                arma::colvec h = AtA * (uvec - vvec);
                                arma::mat hess = AtA * (14 - 105/4 * d + 35/2 *
                                        pow(d,3) - 21/4 * pow(d,5)) -
                                        (105/4/d - 105/2 * d +
                                        105/4 * pow(d,3)) * h * h.t();

                                C(i, j) = arma::as_scalar(d1vec.t() * hess *
                                        d2vec);
                        }
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd1_matern1)]]
NumericMatrix covd1_matern1(NumericMatrix u, NumericMatrix v,
                             NumericMatrix dir1, NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec dvec = D1.row(j).t();

                        double d = arma::as_scalar(sqrt(sum(square(B *
                                                   (uvec - vvec)))));
                        d += 1e-12;

                        arma::colvec h = uvec - vvec;
                        arma::colvec grad = 25 * exp(-5 * d) * AtA * h;

                        C(i, j) = arma::as_scalar(grad.t() * dvec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd2_matern1)]]
NumericMatrix covd2_matern1(NumericMatrix u, NumericMatrix v,
                             NumericMatrix dir1, NumericMatrix dir2,
                             NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);
        arma::mat D2 = arma::mat(dir2.begin(), dir2.nrow(), dir2.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec d1vec = D1.row(i).t();
                        arma::colvec d2vec = D2.row(j).t();

                        double d = arma::as_scalar(sqrt(sum(square(B *
                                                   (uvec - vvec)))));
                        d += 1e-12;

                        arma::colvec h = AtA * (uvec - vvec);
                        arma::mat hess = 25 * exp(-5 * d) * (AtA -
                                5/d * h * h.t());

                        C(i, j) = arma::as_scalar(d1vec.t() * hess * d2vec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd1_matern2)]]
NumericMatrix covd1_matern2(NumericMatrix u, NumericMatrix v,
                            NumericMatrix dir1, NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec dvec = D1.row(j).t();

                        double d = arma::as_scalar(sqrt(sum(square(B *
                                                   (uvec - vvec)))));
                        d += 1e-12;

                        arma::colvec h = uvec - vvec;
                        arma::colvec grad = (12 + 72*d) * exp(-6 * d) * AtA * h;

                        C(i, j) = arma::as_scalar(grad.t() * dvec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd2_matern2)]]
NumericMatrix covd2_matern2(NumericMatrix u, NumericMatrix v,
                            NumericMatrix dir1, NumericMatrix dir2,
                            NumericMatrix A) {
        int N_u = u.nrow();
        int N_v = v.nrow();

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);
        arma::mat D2 = arma::mat(dir2.begin(), dir2.nrow(), dir2.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec d1vec = D1.row(i).t();
                        arma::colvec d2vec = D2.row(j).t();

                        double d = arma::as_scalar(sqrt(sum(square(B *
                                                   (uvec - vvec)))));
                        d += 1e-12;

                        arma::colvec h = AtA * (uvec - vvec);
                        arma::mat hess = (12 + 72*d) * exp(-6 * d) * AtA -
                                432 * exp(-6 * d) * h * h.t();

                        C(i, j) = arma::as_scalar(d1vec.t() * hess * d2vec);
                }
        }

        return(wrap(C));
}


// [[Rcpp::export(.covd1_cauchy)]]
NumericMatrix covd1_cauchy(NumericMatrix u, NumericMatrix v,
                            NumericMatrix dir1, NumericMatrix A,
                            double p) {
        int N_u = u.nrow();
        int N_v = v.nrow();
        double beta = pow(0.05, -1/p) - 1;

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec dvec = D1.row(j).t();

                        double d2 = arma::as_scalar((sum(square(B *
                                                   (uvec - vvec)))));
                        d2 += 1e-12;

                        arma::colvec h = uvec - vvec;
                        arma::colvec grad = 2 * p * beta *
                                pow(1 + beta * d2, -p - 1) * AtA * h;

                        C(i, j) = arma::as_scalar(grad.t() * dvec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.covd2_cauchy)]]
NumericMatrix covd2_cauchy(NumericMatrix u, NumericMatrix v,
                            NumericMatrix dir1, NumericMatrix dir2,
                            NumericMatrix A, double p) {
        int N_u = u.nrow();
        int N_v = v.nrow();
        double beta = pow(0.05, -1/p) - 1;

        arma::mat B = arma::mat(A.begin(), A.nrow(), A.ncol(), false);
        arma::mat AtA = B.t() * B;
        arma::mat U = arma::mat(u.begin(), u.nrow(), u.ncol(), false);
        arma::mat V = arma::mat(v.begin(), v.nrow(), v.ncol(), false);
        arma::mat D1 = arma::mat(dir1.begin(), dir1.nrow(), dir1.ncol(), false);
        arma::mat D2 = arma::mat(dir2.begin(), dir2.nrow(), dir2.ncol(), false);

        NumericMatrix C(N_u, N_v);

        for(int i = 0; i < N_u; i++){
                for(int j = 0; j < N_v; j++){
                        arma::colvec uvec = U.row(i).t();
                        arma::colvec vvec = V.row(j).t();
                        arma::colvec d1vec = D1.row(i).t();
                        arma::colvec d2vec = D2.row(j).t();

                        double d2 = arma::as_scalar((sum(square(B *
                                                    (uvec - vvec)))));
                        d2 += 1e-12;

                        arma::colvec h = AtA * (uvec - vvec);
                        arma::mat hess = 2 * p * beta *
                                pow(1 + beta * d2, -p - 1) * AtA -
                                4 * beta  * p * (p + 1) *
                                pow(1 + beta * d2, -p - 2) * h * h.t();

                        C(i, j) = arma::as_scalar(d1vec.t() * hess * d2vec);
                }
        }

        return(wrap(C));
}

// [[Rcpp::export(.anisotropyC)]]
NumericMatrix anisotropyC(double maxrange, double midrange, double minrange,
                          double azimuth, double dip, double rake){

  static const double pi = 3.14159265;

  // conversion to radians
  azimuth = azimuth * pi / 180;
  dip = dip * pi / 180;
  rake = rake * pi / 180;

  // conversion to mathematical coordinates
  dip = - dip;
  arma::colvec r = arma::colvec(3);
  r(0) = midrange; r(1) = maxrange; r(2) = minrange;

  // rotarion matrices
  arma::mat Rx = arma::mat(3, 3).eye();
  arma::mat Ry = arma::mat(3, 3).eye();
  arma::mat Rz = arma::mat(3, 3).eye();

  Rx(0, 0) = cos(rake);
  Rx(2, 0) = - sin(rake);
  Rx(0, 2) = sin(rake);
  Rx(2, 2) = cos(rake);

  Ry(1, 1) = cos(dip);
  Ry(2, 1) = sin(dip);
  Ry(1, 2) = - sin(dip);
  Ry(2, 2) = cos(dip);

  Rz(0, 0) = cos(azimuth);
  Rz(1, 0) = - sin(azimuth);
  Rz(0, 1) = sin(azimuth);
  Rz(1, 1) = cos(azimuth);

  arma::mat A = arma::diagmat(r);

  A = Rz * Ry * Rx * A;

  return(wrap(A));
}

// [[Rcpp::export(cov_ns)]]
NumericMatrix cov_ns(NumericMatrix x, NumericMatrix y,
                     NumericVector x_sd, NumericVector y_sd,
                     NumericVector x_maxrange, NumericVector y_maxrange,
                     NumericVector x_midrange, NumericVector y_midrange,
                     NumericVector x_minrange, NumericVector y_minrange,
                     NumericVector x_azimuth, NumericVector y_azimuth,
                     NumericVector x_dip, NumericVector y_dip,
                     NumericVector x_rake, NumericVector y_rake,
                     String type, double p = 1){

  NumericMatrix K = NumericMatrix(x.nrow(), y.nrow());

  for(int i = 0; i < x.nrow(); i++){
    NumericMatrix Ax_R = anisotropyC(x_maxrange(i),
                                     x_midrange(i),
                                     x_minrange(i),
                                     x_azimuth(i),
                                     x_dip(i),
                                     x_rake(i));
    arma::mat Ax = arma::mat(Ax_R.begin(), Ax_R.nrow(), Ax_R.ncol(), false);
    Ax = Ax.t() * Ax;

    double Dx = pow(arma::det(Ax), 0.25);

    for(int j = 0; j < y.nrow(); j++){
      NumericMatrix Ay_R = anisotropyC(y_maxrange(j),
                                       y_midrange(j),
                                       y_minrange(j),
                                       y_azimuth(j),
                                       y_dip(j),
                                       y_rake(j));
      arma::mat Ay = arma::mat(Ay_R.begin(), Ay_R.nrow(), Ay_R.ncol(), false);
      Ay = Ay.t() * Ay;

      double Dy = pow(arma::det(Ay), 0.25);

      arma::mat Amix = 0.5 * (Ax + Ay);
      double Dmix = pow(arma::det(Amix), - 0.5);

      // scaled distance
      arma::vec dif = arma::vec(x.row(i) - y.row(j));
      double d = sqrt(arma::as_scalar(dif.t() * arma::solve(Amix, dif)));

      // correlation according to type
      double cor = 0;
      if(type == "gaussian") cor = exp(- 3 * pow(d, 2));
      if(type == "exponential") cor = exp(- 3 * d);
      if(type == "spherical") {
        if(d <= 1) cor = 1 - 1.5 * d + 0.5 * pow(d, 3);
      }
      if(type == "cubic"){
        if(d <= 1) cor =  1 - 7 * pow(d, 2) + 35 / 4 * pow(d, 3) -
          7 / 2 * pow(d, 5) + 3 / 4 * pow(d, 7);
        if(cor < 0) cor = 0; // weird bug gives negative values
      }
      if(type == "matern1") cor = (1 + 5 * d) * exp(- 5 * d);
      if(type == "matern2") cor = (1 + 6 * d + 12 * pow(d, 2)) * exp(- 6 * d);
      if(type == "cauchy") cor = pow(1 + pow(d, 2), - p);

      K(i, j) = Dx * Dy * Dmix * cor * x_sd(i) * y_sd(j);
    }
  }

  return(wrap(K));
}
