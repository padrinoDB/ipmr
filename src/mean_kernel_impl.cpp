#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List mean_kernel_impl(List holder) {

  int n_tot = holder.length();

  NumericMatrix temp = holder[0];

  int n_row = temp.rows();
  int n_col = temp.cols();

  NumericMatrix out(n_row, n_col);

  for(int rw = 0; rw < n_row; rw++) {
    for(int cl = 0; cl < n_col; cl++) {
      for(int i = 0; i < n_tot; i++) {

        NumericMatrix use_mat = holder[i];

        out(rw, cl) += use_mat(rw, cl);

      }
    }
  }

  return List::create(Named("mean_kernel") = out / n_tot);
}
