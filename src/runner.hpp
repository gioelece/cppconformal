#ifndef __RUNNER_HPP
#define __RUNNER_HPP
#include <tuple>
#include <cmath>
#include <optional>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "linear_regr.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd build_grid(Eigen::VectorXd const & ylim, int grid_size);

// [[Rcpp::export]]
List run_linear_conformal(
    Eigen::MatrixXd const & X, Eigen::MatrixXd const & y, Eigen::MatrixXd const & X0,
    int grid_size = 500, double grid_param = 1.25
);
#endif