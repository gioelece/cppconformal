#include "runner.hpp"

MatrixXd build_grid(VectorXd const & ylim, int grid_size) {
    const int d = ylim.size();
    MatrixXd y_grid = MatrixXd::Zero(pow(grid_size, d), d);
    for (int i = 0; i < d; i++) {
        MatrixXd tmp = VectorXd::LinSpaced(grid_size, -ylim(i), ylim(i)) \
            .replicate(pow(grid_size, i), pow(grid_size, d-i-1)).transpose();
        y_grid.col(i) << Map<VectorXd>(tmp.data(), pow(grid_size, d));
    }
    return y_grid;
}

List run_linear_conformal(
    MatrixXd const & X, MatrixXd const & y, MatrixXd const & X0,
    int grid_size, double grid_param
) {
    Eigen::initParallel();
    const int n = X.rows(), n0 = X0.rows(),
              p = X.cols(), d = y.cols();

    // Create a vector containing the grid points that will be tested
    const VectorXd ylim = grid_param * y.array().abs().colwise().maxCoeff();
    const MatrixXd y_grid = build_grid(ylim, grid_size);

    // Create a matrix containing the p-values
    MatrixXd p_values = MatrixXd::Zero(n0, y_grid.rows());
    
    #pragma omp parallel
    {
        // Create a matrix and a vector for the regression model ("y = X \beta + \varepsilon"),
        // adding a row which will contain the point x0 and the corresponding y0 being tested.
        // Its initial value does not really matter, it will be immediately overriden.
        MatrixXd regression_matrix(n + 1, p);
        MatrixXd regression_vector(n + 1, d);
        regression_matrix << X, MatrixXd::Zero(1, p);
        regression_vector << y, MatrixXd::Zero(1, d);

        #pragma omp for collapse(2)
        for(int i = 0; i < n0; i++) {
            for (int j = 0; j < y_grid.rows(); j++) {
                VectorXd y0 = y_grid.row(j);
                regression_matrix.row(n) = X0.row(i);
                regression_vector.row(n) = y0;

                LinearRegression model;
                model.fit(regression_matrix, regression_vector);
                MatrixXd fitted_values = model.predict(regression_matrix);
                ArrayXd residuals = (regression_vector - fitted_values).rowwise().norm().array();

                p_values(i, j) = (residuals > residuals(n)).count() / (n+1.0);
            }
        }
    }

    return List::create(Named("y_grid") = y_grid, 
                        Named("p_values") = p_values);
}
