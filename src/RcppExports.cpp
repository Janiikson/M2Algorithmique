// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_hclust_cpp
List fast_hclust_cpp(NumericMatrix dist_matrix, String method);
RcppExport SEXP _FastHierarchicalClust_fast_hclust_cpp(SEXP dist_matrixSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dist_matrix(dist_matrixSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_hclust_cpp(dist_matrix, method));
    return rcpp_result_gen;
END_RCPP
}
// naive_hclust_cpp
List naive_hclust_cpp(NumericMatrix dist_matrix, String method);
RcppExport SEXP _FastHierarchicalClust_naive_hclust_cpp(SEXP dist_matrixSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dist_matrix(dist_matrixSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(naive_hclust_cpp(dist_matrix, method));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _FastHierarchicalClust_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastHierarchicalClust_fast_hclust_cpp", (DL_FUNC) &_FastHierarchicalClust_fast_hclust_cpp, 2},
    {"_FastHierarchicalClust_naive_hclust_cpp", (DL_FUNC) &_FastHierarchicalClust_naive_hclust_cpp, 2},
    {"_FastHierarchicalClust_rcpp_hello_world", (DL_FUNC) &_FastHierarchicalClust_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastHierarchicalClust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
