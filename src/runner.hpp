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
    Model const & initial_model,
    MatrixXd const & X, MatrixXd const & y, MatrixXd const & X0,
    MatrixXd const & y_grid
);

// [[Rcpp::export]]
List run_linear_conformal_single_grid(
    Eigen::MatrixXd const & X, Eigen::MatrixXd const & y, Eigen::MatrixXd const & X0,
    int grid_size = 500, double grid_param = 1.25
);

Grid create_new_grid_from_pvalues(
    const Grid & old_grid, const VectorXd & p_values, double min_value, int new_grid_size
);

// [[Rcpp::export]]
List run_linear_conformal_multi_grid(
    Eigen::MatrixXd const & X, Eigen::MatrixXd const & y, Eigen::RowVectorXd const & X0,
    Eigen::VectorXd grid_levels, Eigen::VectorXd grid_sizes,
    double initial_grid_param
);
#endif
