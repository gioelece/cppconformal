#include "grid.hpp"

int Grid::get_size() const {
    return pow(grid_points_on_each_axis, d);
}

VectorXd Grid::get_point(int i) const {
    VectorXd point(d);
    auto delta = end_point - start_point;

    for (int axis_idx = 0; axis_idx < d; axis_idx++) {
        int point_idx_on_axis = int(i / pow(grid_points_on_each_axis, axis_idx))
            % grid_points_on_each_axis;
        point(axis_idx) = point_idx_on_axis * delta(axis_idx)
            / (grid_points_on_each_axis - 1)
            + start_point(axis_idx);
    }
    
    return point;
}

MatrixXd Grid::collect() const {
    const int size = get_size();
    MatrixXd y_grid(size, d);
    for (int i = 0; i < size; i++) {
        y_grid.row(i) = get_point(i);
    }
    return y_grid;
}
