#include "grid.hpp"

int Grid::get_size() const {
    return pow(grid_side, d);
}

void Grid::compute_step_increment() {
    step_increment = (end_point - start_point) / (grid_side - 1);
}

VectorXd Grid::get_point(int i) const {
    VectorXd point(d);

    for (int axis_idx = 0; axis_idx < d; axis_idx++) {
        int point_idx_on_axis = int(i / pow(grid_side, axis_idx))
            % grid_side;
        point(axis_idx) = point_idx_on_axis * step_increment(axis_idx)
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
