# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

soft_full <- function(L, lambda) {
    .Call(`_VARDetect_soft_full`, L, lambda)
}

var_break_fit_block_cpp <- function(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index) {
    .Call(`_VARDetect_var_break_fit_block_cpp`, data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index)
}

var_break_fit_block_group_cpp <- function(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index) {
    .Call(`_VARDetect_var_break_fit_block_group_cpp`, data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index)
}

var_break_fit_block_grouprow_cpp <- function(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index) {
    .Call(`_VARDetect_var_break_fit_block_grouprow_cpp`, data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index)
}

var_break_fit_block_groupidx_cpp <- function(data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index) {
    .Call(`_VARDetect_var_break_fit_block_groupidx_cpp`, data, lambda, lambda2, q, max_iteration, tol, initial_phi, blocks, cv_index, group_index)
}

soft_cpp <- function(L, weight, lambda) {
    .Call(`_VARDetect_soft_cpp`, L, weight, lambda)
}

pred_cpp <- function(Y, phi, q, T, k, h) {
    .Call(`_VARDetect_pred_cpp`, Y, phi, q, T, k, h)
}

var_lasso_brk <- function(data, lambda, q, max_iteration, tol) {
    .Call(`_VARDetect_var_lasso_brk`, data, lambda, q, max_iteration, tol)
}

var_lasso_brk_group <- function(data, lambda, q, max_iteration, tol, group_index) {
    .Call(`_VARDetect_var_lasso_brk_group`, data, lambda, q, max_iteration, tol, group_index)
}

var_lasso_brk_group_idx <- function(data, lambda, q, max_iteration, tol, group_index) {
    .Call(`_VARDetect_var_lasso_brk_group_idx`, data, lambda, q, max_iteration, tol, group_index)
}

lambda_warm_up <- function(data, q, blocks, cv_index) {
    .Call(`_VARDetect_lambda_warm_up`, data, q, blocks, cv_index)
}

local_refine <- function(data, q, blocks, nums, lb1, ub2, phi_hat_1, phi_hat_2) {
    .Call(`_VARDetect_local_refine`, data, q, blocks, nums, lb1, ub2, phi_hat_1, phi_hat_2)
}

