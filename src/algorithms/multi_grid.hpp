/*! @file */
#ifndef __ALGORITHMS__MULTI_GRID_HPP
#define __ALGORITHMS__MULTI_GRID_HPP
#include <iostream>
#include <RcppEigen.h>
#include "../grid.hpp"
#include "single_grid.hpp"

/*! Implementation of a multi-grid conformal algorithm.
*   It uses an "inner" single-grid algorithm at each step to recursively select a subgrid.
*/
template<class Model>
class MultiGridAlgorithm : public AlgorithmBase<Model> {
    public:
    /*! Construct a MultiGridAlgorithm instance
        \param grid_levels minimum value of p-values to use at each grid refinement
        \param grid_sides number of points for each side of the grid, for each grid refinement
            (must be one item longer than grid_levels)
        \param initial_grid_param determines the initial size of the grid
        \param inner_algorithm inner algorithm to use (std::unique_ptr)
        \param print_progress print to stdout every run of the inner algorithm.
    */
    MultiGridAlgorithm(
        const VectorXd & _grid_levels, const VectorXd & _grid_sides,
        double _initial_grid_param,
        std::unique_ptr<SingleGridAlgorithm<Model>> _inner_algorithm,
        bool _print_progress = false
    ) : 
        grid_levels(_grid_levels), grid_sides(_grid_sides),
        initial_grid_param(_initial_grid_param),
        print_progress(_print_progress),
        inner_algorithm(std::move(_inner_algorithm))
    {};
        
    /*! Use the p-values from a previous run of the algorithm to identify the areas with p-values greater or equal to min_value,
        and create a new grid;
        \param old_grid grid used to calculate p_values
        \param p_values p-values from the previous run
        \param min_value minimum p-value to consider
        \param new_grid_side number of points for each side of the grid
        \return The new grid object
    */
    static Grid create_new_grid_from_pvalues(
        const Grid & old_grid, const RowVectorXd & p_values, double min_value, int new_grid_side
    );

    /*! Run a conformal algorithm with multi grid refinement,
        computing a confidence region for the covariates corresponding to `Xhat`.

        __Remark__: multi_grid accepts only a single `Xhat`.

        \param model model to use as a base for conformal regression (will be copied at each run)
        \param X matrix of the independent variables
        \param Y matrix of the covariates
        \param Xhat a single point containing the values for the independent variables
        \return An Rcpp list with the following members:
        - `y_grid`: matrix with the coordinates of grid points in the space of the covariates
        - `p_values`: p-values corresponding to those grid points
        - `y_grid_parameters`: vector with the history of grid parameters (start point, end point, grid side)
            for each tried grid
    */
    List run(
        const Model & model,
        const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat
    ) override;

    private:
    VectorXd grid_levels;
    VectorXd grid_sides;
    double initial_grid_param;
    bool print_progress;
    std::unique_ptr<SingleGridAlgorithm<Model>> inner_algorithm;
};

template<class Model>
Grid MultiGridAlgorithm<Model>::create_new_grid_from_pvalues(
    const Grid & old_grid, const RowVectorXd & p_values, double min_value, int new_grid_side
) {
    const ArrayXd old_start = old_grid.get_start_point(), old_end = old_grid.get_end_point(),
            step_increment = old_grid.get_step_increment();
    bool point_found = false;
    ArrayXd start, end;
        
    for (int i = 0; i < old_grid.get_size(); i++) {
        if (p_values(i) >= min_value) {
            ArrayXd coords = old_grid.get_point(i).array();
            if (point_found) {
                start = start.min(coords - step_increment);
                end = end.max(coords + step_increment);
            } else {
                start = coords - step_increment;
                end = coords + step_increment;
            }

            point_found = true;
        }
    }

    if (!point_found) {
        stop("No point over min_value = %d found", min_value);
    }

    return Grid(start.max(old_start), end.min(old_end), new_grid_side);
}


template<class Model>
List MultiGridAlgorithm<Model>::run(
    const Model & model,
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat
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
    grid_parameters.push_back(grid.get_parameters_as_list());

    RowVectorXd p_values;

    int i;
    for (i = 0; i < grid_levels.size(); i++) {
        if (print_progress) {
            std::cout << "Running conformal on grid " << i <<
                " (grid_side = " << grid.get_grid_side() <<
                ", start_point = " << grid.get_start_point().transpose() <<
                ", end_point = " << grid.get_end_point().transpose() <<
                ")" << std::endl;
        }

        p_values = inner_algorithm->run_on_grid(model, X, Y, Xhat, grid);
        grid = create_new_grid_from_pvalues(grid, p_values, grid_levels[i], grid_sides[i+1]);
        grid_parameters.push_back(grid.get_parameters_as_list());
    }

    if (print_progress) {
        std::cout << "Running conformal on grid " << i <<
            " (grid_side = " << grid.get_grid_side() <<
            ", start_point = " << grid.get_start_point().transpose() <<
            ", end_point = " << grid.get_end_point().transpose() <<
            ")" << std::endl;
    }
    p_values = inner_algorithm->run_on_grid(model, X, Y, Xhat, grid);
    
    return List::create(Named("y_grid") = grid.collect(), 
                        Named("y_grid_parameters") = grid_parameters,
                        Named("p_values") = p_values);
}

#endif
