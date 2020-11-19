#include "exports.hpp"
#include "algorithms/multi_grid.hpp"
#include "algorithms/single_grid.hpp"
#include "models/linear_regr.hpp"

List run_linear_conformal_single_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    int grid_side, double grid_param
) {
    LinearRegression model;
    SingleGridAlgorithm<LinearRegression> algorithm(grid_side, grid_param);
    return algorithm.run(model, X, Y, Xhat);
}


List run_ridge_conformal_single_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    double lambda, int grid_side, double grid_param
) {
    RidgeRegression model(lambda);
    SingleGridAlgorithm<RidgeRegression> algorithm(grid_side, grid_param);
    return algorithm.run(model, X, Y, Xhat);
}


List run_linear_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    LinearRegression model;
    SingleGridAlgorithm<LinearRegression> inner_algorithm(grid_sides[0], initial_grid_param);
    MultiGridAlgorithm<LinearRegression> algorithm(grid_levels, grid_sides, initial_grid_param, inner_algorithm, print_progress);
    return algorithm.run(model, X, Y, Xhat);
}


List run_ridge_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat, double lambda,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    RidgeRegression model(lambda);
    SingleGridAlgorithm<RidgeRegression> inner_algorithm(grid_sides[0], initial_grid_param);
    MultiGridAlgorithm<RidgeRegression> algorithm(grid_levels, grid_sides, initial_grid_param, inner_algorithm, print_progress);
    return algorithm.run(model, X, Y, Xhat);
}
