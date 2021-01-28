#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
DataFrame mat_to_df_impl(NumericMatrix mat) {

  int rows = mat.nrow();
  int cols = mat.ncol();

  int out_size = rows * cols;

  NumericVector val(out_size);
  IntegerVector col_ind(out_size);
  IntegerVector row_ind(out_size);

  int it = 0;

  for(int i = 0; i < rows; i++) {

    Rcpp::checkUserInterrupt();

    for(int j = 0; j < cols; j++) {

      val(it) = mat(i, j);
      col_ind(it) = j + 1;
      row_ind(it) = i + 1;
      it = it + 1;

    }

  }

  DataFrame out = DataFrame::create(Named("t") = col_ind,
                                    Named("t_1") = row_ind,
                                    Named("value") = val);

  return out;

}
