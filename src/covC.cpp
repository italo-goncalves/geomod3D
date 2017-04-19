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
