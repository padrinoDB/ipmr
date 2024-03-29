// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mat_to_df_impl
DataFrame mat_to_df_impl(NumericMatrix mat);
RcppExport SEXP _ipmr_mat_to_df_impl(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_to_df_impl(mat));
    return rcpp_result_gen;
END_RCPP
}
// mean_kernel_impl
List mean_kernel_impl(List holder);
RcppExport SEXP _ipmr_mean_kernel_impl(SEXP holderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type holder(holderSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_kernel_impl(holder));
    return rcpp_result_gen;
END_RCPP
}
// update_pop_state
List update_pop_state(List holder, List current, int iter);
RcppExport SEXP _ipmr_update_pop_state(SEXP holderSEXP, SEXP currentSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type holder(holderSEXP);
    Rcpp::traits::input_parameter< List >::type current(currentSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(update_pop_state(holder, current, iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ipmr_mat_to_df_impl", (DL_FUNC) &_ipmr_mat_to_df_impl, 1},
    {"_ipmr_mean_kernel_impl", (DL_FUNC) &_ipmr_mean_kernel_impl, 1},
    {"_ipmr_update_pop_state", (DL_FUNC) &_ipmr_update_pop_state, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ipmr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
