#include "runner.hpp"
#include "Eigen/src/Core/VectorBlock.h"

MatrixXd run_linear_conformal(
    MatrixXd const & X, VectorXd const & y, MatrixXd const & X0,
    int grid_size, double grid_param
) {
    Eigen::initParallel();
    int n = X.rows(), n0 = X0.rows(),
        p = X.cols();

    // Create a vector containing the grid points that will be tested
    double ylim = grid_param * y.array().abs().maxCoeff();
    VectorXd y_grid = VectorXd::LinSpaced(grid_size, -ylim, ylim);

    // Create a matrix containing the p-values
    MatrixXd p_values = MatrixXd::Zero(X0.rows(), grid_size);
    
    #pragma omp parallel
    {
        // Create a matrix and a vector for the regression model ("y = X \beta + \varepsilon"),
        // adding a row which will contain the point x0 and the corresponding y0 being tested.
        // Its initial value does not really matter, it will be immediately overriden.
        MatrixXd regression_matrix(n + 1, p);
        VectorXd regression_vector(n + 1);
        regression_matrix << X, X0.row(0);
        regression_vector << y, 0;

        #pragma omp for collapse(2)
        for(int i = 0; i < n0; i++) {
            for (int j = 0; j < grid_size; j++) {
                double y0 = y_grid(j);
                #pragma omp critical
                std::cout << i << " " << j << "\n";
                regression_matrix.row(n) = X0.row(i);
                regression_vector(n) = y0;

                LinearRegression model;
                model.fit(regression_matrix, regression_vector);
                VectorXd fitted_values = model.predict(regression_matrix);
                ArrayXd residuals = (regression_vector - fitted_values).array().abs();

                std::cout << i << " " << j << "\n";
                p_values(i, j) = (residuals > residuals(n)).count() / (n+1.0);
            }
        }
    }
    return p_values;
}
