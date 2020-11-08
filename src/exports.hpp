/*! @file */
#ifndef __EXPORTS_HHP
#define __EXPORTS_HHP
#include "algorithms/multi_grid.hpp"
#include "algorithms/single_grid.hpp"
#include "models/linear_regr.hpp"

// Remark: Rcpp does not work if the Eigen namespace is omitted from exported definitions.

// [[Rcpp::export]]
/*! Run a conformal algorithm with a simple grid and a linear regression model.
    For details, see @ref SingleGridAlgorithm::run.
*/
List run_linear_conformal_single_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat,
    int grid_side = 500, double grid_param = 1.25
);

// [[Rcpp::export]]
/*! Run a conformal algorithm with a simple grid and a ridge regression model.
    For details, see @ref SingleGridAlgorithm::run.

    \param lambda lambda parameter for the ridge regression
*/
List run_ridge_conformal_single_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat,
    double lambda, int grid_side = 500, double grid_param = 1.25
);


/*! Run a conformal algorithm with automatic multi grid refinement and linear regression model.
    See @ref MultiGridAlgorithm::run for details.
*/
// [[Rcpp::export]]
List run_linear_conformal_multi_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param,
    bool print_progress = false
);

/*! Run a conformal algorithm with automatic multi grid refinement and ridge regression model.
    See @ref MultiGridAlgorithm::run for details.

    \param lambda lambda parameter for the ridge regression
*/
// [[Rcpp::export]]
List run_ridge_conformal_multi_grid(
    const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat, double lambda,
    const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param,
    bool print_progress = false
);

#endif
