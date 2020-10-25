#include "runner.hpp"

template<class Model>
MatrixXd run_conformal_on_grid(
    Model const & initial_model,
    MatrixXd const & X, MatrixXd const & y, MatrixXd const & X0,
    Grid const & grid
) {
    const int n = X.rows(), n0 = X0.rows(),
              p = X.cols(), d = y.cols(),
              num_threads = omp_get_max_threads();
    // Create a matrix containing the p-values
    MatrixXd p_values = MatrixXd::Zero(n0, grid.get_size());
    // Create a matrix and a vector for the regression model ("y = X \beta + \varepsilon"),
    // adding (for each thread) a row which will contain the point x0 and the corresponding y0 being tested.
    // Its initial value does not really matter, it will be immediately overriden.
    MatrixXd regression_matrix(n + num_threads, p);
    MatrixXd regression_vector(n + num_threads, d);
    regression_matrix << X, MatrixXd::Zero(num_threads, p);
    regression_vector << y, MatrixXd::Zero(num_threads, d);

    #pragma omp parallel
    {
        int this_thread = omp_get_thread_num();
        SparseMatrix<double> reduction_matrix(n + 1, n + num_threads);
        reduction_matrix.reserve(VectorXi::Ones(n + num_threads));
        for (int i = 0; i < n; i++) {
            reduction_matrix.insert(i, i) = 1;
        }
        reduction_matrix.insert(n, n + this_thread) = 1;

        #pragma omp for collapse(2)
        for(int i = 0; i < n0; i++) {
            for (int j = 0; j < grid.get_size(); j++) {
                VectorXd y0 = grid.get_point(j);
                regression_matrix.row(n + this_thread) = X0.row(i);
                regression_vector.row(n + this_thread) = y0;

                Model model(initial_model);
                auto thread_regression_matrix = reduction_matrix * regression_matrix;
                auto thread_regression_vector = reduction_matrix * regression_vector;
                
                model.fit(thread_regression_matrix, thread_regression_vector);
                MatrixXd fitted_values = model.predict(thread_regression_matrix);
                ArrayXd residuals = (thread_regression_vector - fitted_values).rowwise().norm().array();

                p_values(i, j) = (residuals > residuals(n)).count() / (n+1.0);
            }
        }
    }

    return p_values;
}


List run_linear_conformal(
    MatrixXd const & X, MatrixXd const & y, MatrixXd const & X0,
    int grid_size, double grid_param
) {
    // Create a vector containing the grid points that will be tested
    const VectorXd ylim = grid_param * y.array().abs().colwise().maxCoeff();
    const Grid grid(-ylim, ylim, grid_size);

    LinearRegression model;
    MatrixXd p_values = run_conformal_on_grid(model, X, y, X0, grid);

    return List::create(Named("y_grid") = grid.collect(), 
                        Named("p_values") = p_values);
}
