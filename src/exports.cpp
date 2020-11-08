#include "exports.hpp"
#include "algorithms/multi_grid.hpp"

List run_linear_conformal_single_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    int grid_side, double grid_param
) {
    LinearRegression model;
    SingleGridAlgorithm algorithm;
    return algorithm.run(model, X, Y, Xhat, grid_side, grid_param);
}


List run_ridge_conformal_single_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    double lambda, int grid_side, double grid_param
) {
    RidgeRegression model(lambda);
    SingleGridAlgorithm algorithm;
    return algorithm.run(model, X, Y, Xhat, grid_side, grid_param);
}


List run_linear_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    LinearRegression model;
    MultiGridAlgorithm algorithm;
    return algorithm.run(model, X, Y, Xhat, grid_levels, grid_sides, initial_grid_param, print_progress);
}


List run_ridge_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat, double lambda,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    RidgeRegression model(lambda);
    MultiGridAlgorithm algorithm;
    return algorithm.run(model, X, Y, Xhat, grid_levels, grid_sides, initial_grid_param, print_progress);
}
