#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::VectorXd right_mult_impl(Eigen::Map<Eigen::MatrixXd> k,
                                Eigen::Map<Eigen::VectorXd> v) {

  Eigen::VectorXd out = k * v;

  return(out);

}

// [[Rcpp::export]]

Eigen::VectorXd left_mult_impl(Eigen::Map<Eigen::MatrixXd> k,
                               Eigen::Map<Eigen::VectorXd> v) {

  Eigen::VectorXd out = k.transpose() * v;

  return(out);


}
