#ifndef __RUNNER_HPP
#define __RUNNER_HPP
#include <tuple>
#include <cmath>
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "grid.hpp"
#include "linear_regr.hpp"

using namespace Rcpp;
using namespace Eigen;

// Remark: Rcpp does not work if the Eigen namespace is omitted from exported definitions.

template<class Model>
MatrixXd run_conformal_on_grid(
    const Model & initial_model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    const Grid & grid
);

template<class Model>
List run_conformal_single_grid(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    int grid_side, double grid_param
);

// [[Rcpp::export]]
List run_linear_conformal_single_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0,
    int grid_side = 500, double grid_param = 1.25
);

// [[Rcpp::export]]
List run_ridge_conformal_single_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0,
    double lambda, int grid_side = 500, double grid_param = 1.25
);

Grid create_new_grid_from_pvalues(
    const Grid & old_grid, const RowVectorXd & p_values, double min_value, int new_grid_side
);

// Remark: multi_grid accepts only a single X0
template<class Model>
List run_conformal_multi_grid(
    const Model & model,
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::RowVectorXd & X0,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param
);

// [[Rcpp::export]]
List run_linear_conformal_multi_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param
);

// [[Rcpp::export]]
List run_ridge_conformal_multi_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0, double lambda,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param
);
#endif
