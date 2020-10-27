#include "runner.hpp"
#include "linear_regr.hpp"

template<class Model>
MatrixXd run_conformal_on_grid(
    const Model & initial_model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    const Grid & grid
) {
    if (X.cols() != X0.cols()) {
        stop("X.cols() != X0.cols(), but they must be equal (to p)");
    }
    if (X.rows() != Y.rows()) {
        stop("X.rows() != y.rows(), but they must be equal (to n)");
    }

    const int n = X.rows(), n0 = X0.rows(),
              p = X.cols(), d = Y.cols(),
              num_threads = omp_get_max_threads();
    // Create a matrix containing the p-values
    MatrixXd p_values = MatrixXd::Zero(n0, grid.get_size());
    // Create a matrix and a vector for the regression model ("y = X \beta + \varepsilon"),
    // adding (for each thread) a row which will contain the point x0 and the corresponding y0 being tested.
    // Its initial value does not really matter, it will be immediately overriden.
    MatrixXd regression_matrix(n + num_threads, p);
    MatrixXd regression_vector(n + num_threads, d);
    regression_matrix << X, MatrixXd::Zero(num_threads, p);
    regression_vector << Y, MatrixXd::Zero(num_threads, d);

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


template<class Model>
List run_conformal_single_grid(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    int grid_side, double grid_param
) {
    const VectorXd ylim = grid_param * Y.array().abs().colwise().maxCoeff();
    const Grid grid(-ylim, ylim, grid_side);
    MatrixXd p_values = run_conformal_on_grid(model, X, Y, X0, grid);

    return List::create(Named("y_grid") = grid.collect(), 
                        Named("p_values") = p_values);
}


List run_linear_conformal_single_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    int grid_side, double grid_param
) {
    LinearRegression model;
    return run_conformal_single_grid(model, X, Y, X0, grid_side, grid_param);
}


List run_ridge_conformal_single_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    double lambda, int grid_side, double grid_param
) {
    RidgeRegression model(lambda);
    return run_conformal_single_grid(model, X, Y, X0, grid_side, grid_param);
}


Grid create_new_grid_from_pvalues(
    const Grid & old_grid, const RowVectorXd & p_values, double min_value, int new_grid_side
) {
    // TODO: check that there is at least a point with p >= min_value
    ArrayXd start = old_grid.get_end_point(), end = old_grid.get_start_point(),
            step_increment = old_grid.get_step_increment();
        
    for (int i = 0; i < old_grid.get_size(); i++) {
        if (p_values(i) >= min_value) {
            start = start.min(old_grid.get_point(i).array() - step_increment);
            end = end.max(old_grid.get_point(i).array() + step_increment);
        }
    }

    if ((start == old_grid.get_end_point().array()).all() && (end == old_grid.get_start_point().array()).all()) {
        stop("No point over min_value = %d found", min_value);
    }

    return Grid(start, end, new_grid_side);
}


template<class Model>
List run_conformal_multi_grid(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param
) {
    if (X0.rows() > 1) {
        stop("You must pass a single X0 point to multi_grid functions");
    }

    const VectorXd initial_ylim = initial_grid_param * Y.array().abs().colwise().maxCoeff();
    Grid grid(-initial_ylim, initial_ylim, grid_sides[0]);

    RowVectorXd p_values;

    for (int i = 0; i < grid_levels.size(); i++) {
        p_values = run_conformal_on_grid(model, X, Y, X0, grid);
        grid = create_new_grid_from_pvalues(grid, p_values, grid_levels[i], grid_sides[i+1]);
    }

    return List::create(Named("y_grid") = grid.collect(), 
                        Named("p_values") = p_values);
}


List run_linear_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param
) {
    LinearRegression model;
    return run_conformal_multi_grid(model, X, Y, X0, grid_levels, grid_sides, initial_grid_param);
}


List run_ridge_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & X0, double lambda,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param
) {
    RidgeRegression model(lambda);
    return run_conformal_multi_grid(model, X, Y, X0, grid_levels, grid_sides, initial_grid_param);
}
