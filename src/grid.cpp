#include "grid.hpp"

int Grid::get_size() const {
    return pow(grid_side, d);
}

void Grid::compute_step_increment() {
    step_increment = (end_point - start_point) / (grid_side - 1);
}

VectorXd Grid::get_point(int point_idx) const {
    VectorXd point(d);

    for (int i = 0; i < d; i++) {
      int point_idx_on_side = int(point_idx / pow(grid_side, i)) % grid_side;
      point(i) = point_idx_on_side * step_increment(i) + start_point(i);
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
