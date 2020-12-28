#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List update_pop_state(List holder, List current, int iter){

  for(int i = 0; i < holder.size(); i++) {

    NumericVector insert(current[i]);
    NumericMatrix update_mat = holder[i];
    NumericMatrix::Column update_col = update_mat(_, iter);

    update_col = insert;

  }

  return holder;

}
