#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List mean_kernel_impl(List holder) {

  int n_tot = holder.length();

  NumericMatrix temp = holder[0];

  int n_row = temp.rows();
  int n_col = temp.cols();

  NumericMatrix out(n_row, n_col);


  for(int i = 0; i < n_tot; i++) {

    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    NumericMatrix use_mat = holder[i];

    out += use_mat;

  }


  return List::create(Named("mean_kernel") = out / n_tot);
}
