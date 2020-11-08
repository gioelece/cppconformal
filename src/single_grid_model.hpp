/*! @file */
#ifndef __SINGLE_GRID_HPP
#define __SINGLE_GRID_HPP
#include <cmath>
#include <omp.h>
#include <RcppEigen.h>
#include "grid.hpp"
#include "linear_regr.hpp"

using namespace Rcpp;
using namespace Eigen;

// Remark: Rcpp does not work if the Eigen namespace is omitted from exported definitions.

/*! Run a conformal algorithm on a @ref Grid instance.

    \param model model to use as a base for conformal regression (will be copied at each run)
    \param X matrix of the independent variables
    \param Y matrix of the covariates
    \param Xhat a matrix containing multiple points to use as values for the independent variables
    \param grid grid instance
    \return An Rcpp list with the following members:
    - `y_grid`: matrix with the coordinates of grid points in the space of the covariates
    - `p_values`: p-values corresponding to those grid points
*/
template<class Model>
MatrixXd run_conformal_on_grid(
    const Model & initial_model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    const Grid & grid
);

/*! Run a conformal algorithm with a simple grid,
    computing a confidence region for the covariates corresponding to `Xhat`.

    \param model model to use as a base for conformal regression (will be copied at each run)
    \param X matrix of the independent variables
    \param Y matrix of the covariates
    \param Xhat a matrix containing multiple points to use as values for the independent variables
    \param grid_side number of points for each side of the grid
    \param initial_grid_param determines the initial size of the grid
    \return An Rcpp list with the following members:
    - `y_grid`: matrix with the coordinates of grid points in the space of the covariates
    - `p_values`: p-values corresponding to those grid points
*/
template<class Model>
List run_conformal_single_grid(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    int grid_side, double grid_param
);

// [[Rcpp::export]]
/*! Run a conformal algorithm with a simple grid and a linear regression model.
    For details, see @ref run_conformal_single_grid.
*/
List run_linear_conformal_single_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat,
    int grid_side = 500, double grid_param = 1.25
);

// [[Rcpp::export]]
/*! Run a conformal algorithm with a simple grid and a ridge regression model.
    For details, see @ref run_conformal_single_grid.

    \param lambda lambda parameter for the ridge regression
*/
List run_ridge_conformal_single_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat,
    double lambda, int grid_side = 500, double grid_param = 1.25
);
#endif
