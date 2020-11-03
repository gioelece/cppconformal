/*! @file */
#ifndef __MULTI_GRID_MODEL_HPP
#define __MULTI_GRID_MODEL_HPP
#include <RcppEigen.h>
#include "grid.hpp"
#include "linear_regr.hpp"
#include "single_grid_model.hpp"

using namespace Rcpp;
using namespace Eigen;

// Remark: Rcpp does not work if the Eigen namespace is omitted from exported definitions.

/*! Use the p-values from a previous run of the algorithm to identify the areas with p-values greater or equal to min_value,
    and create a new grid;
    \param old_grid grid used to calculate p_values
    \param p_values p-values from the previous run
    \param min_value minimum p-value to consider
    \param new_grid_side number of points for each side of the grid
    \return The new grid object
*/
Grid create_new_grid_from_pvalues(
    const Grid & old_grid, const RowVectorXd & p_values, double min_value, int new_grid_side
);

/*! Run a conformal algorithm with multi grid refinement,
    computing a confidence region for the covariates corresponding to `X0`.

    __Remark__: multi_grid accepts only a single `X0`

    \param model model to use as a base for conformal regression (will be copied at each run)
    \param X matrix of the independent variables
    \param Y matrix of the covariates
    \param X0 a single point containing the values for the independent variables
    \param grid_levels minimum value of p-values to use at each grid refinement
    \param grid_sides number of points for each side of the grid, for each grid refinement
    \param initial_grid_param determines the initial size of the grid
    \return An Rcpp list with the following members:
    - `y_grid`: matrix with the coordinates of grid points in the space of the covariates
    - `p_values`: p-values corresponding to those grid points
*/
template<class Model>
List run_conformal_multi_grid(
    const Model & model,
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param
);

/*! Run a conformal algorithm with automatic multi grid refinement and linear regression model.
    See @ref run_conformal_multi_grid for details.
*/
// [[Rcpp::export]]
List run_linear_conformal_multi_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param
);

/*! Run a conformal algorithm with automatic multi grid refinement and ridge regression model.
    See @ref run_conformal_multi_grid for details.

    \param lambda lambda parameter for the ridge regression
*/
// [[Rcpp::export]]
List run_ridge_conformal_multi_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & X0, double lambda,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param
);
#endif
