/*! @file */
#ifndef __ALGORITHMS__SINGLE_GRID_HPP
#define __ALGORITHMS__SINGLE_GRID_HPP
#include <cmath>
#include <omp.h>
#include <RcppEigen.h>
#include "../grid.hpp"

using Rcpp::stop;
using Rcpp::Named;

class SingleGridAlgorithm {
    public:
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
    MatrixXd run_on_grid(
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
    List run(
        const Model & model,
        const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
        int grid_side, double grid_param
    );
};


template<class Model>
MatrixXd SingleGridAlgorithm::run_on_grid(
    const Model & initial_model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    const Grid & grid
) {
    if (X.cols() != Xhat.cols()) {
        stop("X.cols() != Xhat.cols(), but they must be equal (to p)");
    }
    if (X.rows() != Y.rows()) {
        stop("X.rows() != y.rows(), but they must be equal (to n)");
    }

    const int n = X.rows(), n0 = Xhat.rows(),
              p = X.cols(), d = Y.cols(),
              num_threads = omp_get_max_threads();
    // Create a matrix containing the p-values
    MatrixXd p_values = MatrixXd::Zero(n0, grid.get_size());
    // Create a matrix and a vector for the regression model ("y = X \beta + \varepsilon"),
    // adding (for each thread) a row which will contain the point Xhat and the corresponding y0 being tested.
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
                regression_matrix.row(n + this_thread) = Xhat.row(i);
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
List SingleGridAlgorithm::run(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    int grid_side, double grid_param
) {
    const VectorXd ylim = grid_param * Y.array().abs().colwise().maxCoeff();
    const Grid grid(-ylim, ylim, grid_side);
    MatrixXd p_values = run_on_grid(model, X, Y, Xhat, grid);

    return List::create(Named("y_grid") = grid.collect(), 
                        Named("p_values") = p_values);
}

#endif
