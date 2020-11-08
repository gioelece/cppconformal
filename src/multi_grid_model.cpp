#include "multi_grid_model.hpp"

Grid create_new_grid_from_pvalues(
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


template<class Model>
List run_conformal_multi_grid(
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

    RowVectorXd p_values;

    int i;
    for (i = 0; i < grid_levels.size(); i++) {
        if (print_progress) {
            std::cout << "Running conformal on grid " << i <<
                " (grid_side = " << grid.get_grid_side() << ")" << std::endl;
        }

        p_values = run_conformal_on_grid(model, X, Y, Xhat, grid);
        grid = create_new_grid_from_pvalues(grid, p_values, grid_levels[i], grid_sides[i+1]);
        grid_parameters.push_back(grid.get_parameters_as_list());
    }

    std::cout << "Running conformal on grid " << i <<
                " (grid_side = " << grid.get_grid_side() << ")" << std::endl;
    p_values = run_conformal_on_grid(model, X, Y, Xhat, grid);
    
    return List::create(Named("y_grid") = grid.collect(), 
                        Named("y_grid_parameters") = grid_parameters,
                        Named("p_values") = p_values);
}


List run_linear_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    LinearRegression model;
    return run_conformal_multi_grid(model, X, Y, Xhat, grid_levels, grid_sides, initial_grid_param, print_progress);
}


List run_ridge_conformal_multi_grid(
    const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat, double lambda,
    const VectorXd & grid_levels, const VectorXd & grid_sides, double initial_grid_param,
    bool print_progress
) {
    RidgeRegression model(lambda);
    return run_conformal_multi_grid(model, X, Y, Xhat, grid_levels, grid_sides, initial_grid_param, print_progress);
}
