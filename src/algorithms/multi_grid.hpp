/*! @file */
#ifndef __ALGORITHMS__MULTI_GRID_HPP
#define __ALGORITHMS__MULTI_GRID_HPP
#include <iostream>
#include <RcppEigen.h>
#include "../grid.hpp"
#include "single_grid.hpp"

/*! Implementation of a multi-grid conformal algorithm
*/
class MultiGridAlgorithm {
    public:
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
        computing a confidence region for the covariates corresponding to `Xhat`.

        __Remark__: multi_grid accepts only a single `Xhat`.

        \param InnerAlgorithm class to use as inner, single-grid algorithm
        \param model model to use as a base for conformal regression (will be copied at each run)
        \param X matrix of the independent variables
        \param Y matrix of the covariates
        \param Xhat a single point containing the values for the independent variables
        \param grid_levels minimum value of p-values to use at each grid refinement
        \param grid_sides number of points for each side of the grid, for each grid refinement
            (must be one item longer than grid_levels)
        \param initial_grid_param determines the initial size of the grid
        \return An Rcpp list with the following members:
        - `y_grid`: matrix with the coordinates of grid points in the space of the covariates
        - `p_values`: p-values corresponding to those grid points
        - `y_grid_parameters`: vector with the history of grid parameters (start point, end point, grid side)
            for each tried grid
    */
    template<class Model, class InnerAlgorithm = SingleGridAlgorithm>
    List run(
        const Model & model,
        const Eigen::MatrixXd & X, const Eigen::MatrixXd & Y, const Eigen::MatrixXd & Xhat,
        const Eigen::VectorXd & grid_levels, const Eigen::VectorXd & grid_sides, double initial_grid_param,
        bool print_progress = false
    );
};


Grid MultiGridAlgorithm::create_new_grid_from_pvalues(
    const Grid & old_grid, const RowVectorXd & p_values, double min_value, int new_grid_side
) {
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


template<class Model, class InnerAlgorithm>
List MultiGridAlgorithm::run(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    if (Xhat.rows() > 1) {
        stop("You must pass a single Xhat point to multi_grid functions");
    }

    if (grid_levels.size() + 1 != grid_sides.size()) {
        stop("grid_sides must be one item longer than grid_levels");
    }

    std::vector<List> grid_parameters;
    const VectorXd initial_ylim = initial_grid_param * Y.array().abs().colwise().maxCoeff();
    Grid grid(-initial_ylim, initial_ylim, grid_sides[0]);

    InnerAlgorithm inner_algorithm;
    RowVectorXd p_values;

    int i;
    for (i = 0; i < grid_levels.size(); i++) {
        if (print_progress) {
            std::cout << "Running conformal on grid " << i <<
                " (grid_side = " << grid.get_grid_side() << ")" << std::endl;
        }

        p_values = inner_algorithm.run_on_grid(model, X, Y, Xhat, grid);
        grid = create_new_grid_from_pvalues(grid, p_values, grid_levels[i], grid_sides[i+1]);
        grid_parameters.push_back(grid.get_parameters_as_list());
    }

    std::cout << "Running conformal on grid " << i <<
                " (grid_side = " << grid.get_grid_side() << ")" << std::endl;
    p_values = inner_algorithm.run_on_grid(model, X, Y, Xhat, grid);
    
    return List::create(Named("y_grid") = grid.collect(), 
                        Named("y_grid_parameters") = grid_parameters,
                        Named("p_values") = p_values);
}

#endif
