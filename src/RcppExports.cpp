// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// soft_full
arma::mat soft_full(arma::mat L, double lambda);
RcppExport SEXP _VARDetect_soft_full(SEXP LSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_full(L, lambda));
    return rcpp_result_gen;
END_RCPP
}
// var_break_fit_block_cpp
List var_break_fit_block_cpp(NumericMatrix data, double lambda, double lambda2, int q, int max_iteration, double tol, NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index);
RcppExport SEXP _VARDetect_var_break_fit_block_cpp(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lambda2SEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP, SEXP initial_phiSEXP, SEXP blocksSEXP, SEXP cv_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initial_phi(initial_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv_index(cv_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(var_break_fit_block_cpp(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index));
    return rcpp_result_gen;
END_RCPP
}
// var_break_fit_block_group_cpp
List var_break_fit_block_group_cpp(NumericMatrix data, double lambda, double lambda2, int q, int max_iteration, double tol, NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index, List group_index);
RcppExport SEXP _VARDetect_var_break_fit_block_group_cpp(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lambda2SEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP, SEXP initial_phiSEXP, SEXP blocksSEXP, SEXP cv_indexSEXP, SEXP group_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initial_phi(initial_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv_index(cv_indexSEXP);
    Rcpp::traits::input_parameter< List >::type group_index(group_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(var_break_fit_block_group_cpp(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index));
    return rcpp_result_gen;
END_RCPP
}
// var_break_fit_block_grouprow_cpp
List var_break_fit_block_grouprow_cpp(NumericMatrix data, double lambda, double lambda2, int q, int max_iteration, double tol, NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index, List group_index);
RcppExport SEXP _VARDetect_var_break_fit_block_grouprow_cpp(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lambda2SEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP, SEXP initial_phiSEXP, SEXP blocksSEXP, SEXP cv_indexSEXP, SEXP group_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initial_phi(initial_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv_index(cv_indexSEXP);
    Rcpp::traits::input_parameter< List >::type group_index(group_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(var_break_fit_block_grouprow_cpp(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index));
    return rcpp_result_gen;
END_RCPP
}
// var_break_fit_block_groupidx_cpp
List var_break_fit_block_groupidx_cpp(NumericMatrix data, double lambda, double lambda2, int q, int max_iteration, double tol, NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index, List group_index);
RcppExport SEXP _VARDetect_var_break_fit_block_groupidx_cpp(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lambda2SEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP, SEXP initial_phiSEXP, SEXP blocksSEXP, SEXP cv_indexSEXP, SEXP group_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initial_phi(initial_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv_index(cv_indexSEXP);
    Rcpp::traits::input_parameter< List >::type group_index(group_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(var_break_fit_block_groupidx_cpp(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index));
    return rcpp_result_gen;
END_RCPP
}
// soft_cpp
arma::mat soft_cpp(arma::mat L, arma::vec weight, double lambda);
RcppExport SEXP _VARDetect_soft_cpp(SEXP LSEXP, SEXP weightSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_cpp(L, weight, lambda));
    return rcpp_result_gen;
END_RCPP
}
// pred_cpp
arma::mat pred_cpp(arma::mat Y, arma::mat phi, int q, int T, int k, int h);
RcppExport SEXP _VARDetect_pred_cpp(SEXP YSEXP, SEXP phiSEXP, SEXP qSEXP, SEXP TSEXP, SEXP kSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_cpp(Y, phi, q, T, k, h));
    return rcpp_result_gen;
END_RCPP
}
// var_lasso_brk
List var_lasso_brk(NumericMatrix data, NumericVector lambda, int q, int max_iteration, double tol);
RcppExport SEXP _VARDetect_var_lasso_brk(SEXP dataSEXP, SEXP lambdaSEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lasso_brk(data, lambda, q, max_iteration, tol));
    return rcpp_result_gen;
END_RCPP
}
// var_lasso_brk_group
List var_lasso_brk_group(NumericMatrix data, NumericVector lambda, int q, int max_iteration, double tol, List group_index);
RcppExport SEXP _VARDetect_var_lasso_brk_group(SEXP dataSEXP, SEXP lambdaSEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP, SEXP group_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< List >::type group_index(group_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lasso_brk_group(data, lambda, q, max_iteration, tol, group_index));
    return rcpp_result_gen;
END_RCPP
}
// var_lasso_brk_group_idx
List var_lasso_brk_group_idx(NumericMatrix data, NumericVector lambda, int q, int max_iteration, double tol, List group_index);
RcppExport SEXP _VARDetect_var_lasso_brk_group_idx(SEXP dataSEXP, SEXP lambdaSEXP, SEXP qSEXP, SEXP max_iterationSEXP, SEXP tolSEXP, SEXP group_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iteration(max_iterationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< List >::type group_index(group_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(var_lasso_brk_group_idx(data, lambda, q, max_iteration, tol, group_index));
    return rcpp_result_gen;
END_RCPP
}
// lambda_warm_up
List lambda_warm_up(NumericMatrix data, int q, NumericVector blocks, NumericVector cv_index);
RcppExport SEXP _VARDetect_lambda_warm_up(SEXP dataSEXP, SEXP qSEXP, SEXP blocksSEXP, SEXP cv_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv_index(cv_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_warm_up(data, q, blocks, cv_index));
    return rcpp_result_gen;
END_RCPP
}
// local_refine
List local_refine(NumericMatrix data, int q, NumericVector blocks, NumericVector nums, int lb1, int ub2, NumericMatrix phi_hat_1, NumericMatrix phi_hat_2);
RcppExport SEXP _VARDetect_local_refine(SEXP dataSEXP, SEXP qSEXP, SEXP blocksSEXP, SEXP numsSEXP, SEXP lb1SEXP, SEXP ub2SEXP, SEXP phi_hat_1SEXP, SEXP phi_hat_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nums(numsSEXP);
    Rcpp::traits::input_parameter< int >::type lb1(lb1SEXP);
    Rcpp::traits::input_parameter< int >::type ub2(ub2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_hat_1(phi_hat_1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi_hat_2(phi_hat_2SEXP);
    rcpp_result_gen = Rcpp::wrap(local_refine(data, q, blocks, nums, lb1, ub2, phi_hat_1, phi_hat_2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VARDetect_soft_full", (DL_FUNC) &_VARDetect_soft_full, 2},
    {"_VARDetect_var_break_fit_block_cpp", (DL_FUNC) &_VARDetect_var_break_fit_block_cpp, 9},
    {"_VARDetect_var_break_fit_block_group_cpp", (DL_FUNC) &_VARDetect_var_break_fit_block_group_cpp, 10},
    {"_VARDetect_var_break_fit_block_grouprow_cpp", (DL_FUNC) &_VARDetect_var_break_fit_block_grouprow_cpp, 10},
    {"_VARDetect_var_break_fit_block_groupidx_cpp", (DL_FUNC) &_VARDetect_var_break_fit_block_groupidx_cpp, 10},
    {"_VARDetect_soft_cpp", (DL_FUNC) &_VARDetect_soft_cpp, 3},
    {"_VARDetect_pred_cpp", (DL_FUNC) &_VARDetect_pred_cpp, 6},
    {"_VARDetect_var_lasso_brk", (DL_FUNC) &_VARDetect_var_lasso_brk, 5},
    {"_VARDetect_var_lasso_brk_group", (DL_FUNC) &_VARDetect_var_lasso_brk_group, 6},
    {"_VARDetect_var_lasso_brk_group_idx", (DL_FUNC) &_VARDetect_var_lasso_brk_group_idx, 6},
    {"_VARDetect_lambda_warm_up", (DL_FUNC) &_VARDetect_lambda_warm_up, 4},
    {"_VARDetect_local_refine", (DL_FUNC) &_VARDetect_local_refine, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_VARDetect(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
