#ifndef __RUNNER_HPP
#define __RUNNER_HPP
#include <tuple>
#include <cmath>
#include <optional>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "linear_regr.hpp"

using namespace Rcpp;
using namespace Eigen;

// Remark: Rcpp does not work if the Eigen namespace is removed from exported definitions.

// [[Rcpp::export]]
Eigen::MatrixXd build_grid(Eigen::VectorXd const & ylim, int grid_size);

template<class Model>
MatrixXd run_conformal_on_grid(
    Model const & initial_model,
    MatrixXd const & X, MatrixXd const & y, MatrixXd const & X0,
    MatrixXd const & y_grid
);

// [[Rcpp::export]]
List run_linear_conformal(
    Eigen::MatrixXd const & X, Eigen::MatrixXd const & y, Eigen::MatrixXd const & X0,
    int grid_size = 500, double grid_param = 1.25
);
#endif
