#ifndef __GRID_HPP
#define __GRID_HPP
#include <vector>
#include <Eigen/Core>

using namespace Eigen;

class Grid {
    public:
    Grid(const VectorXd & s, const VectorXd & e, int g) :
        start_point(s), end_point(e), grid_points_on_each_axis(g),
        d(s.size()) {};
    int get_size() const;
    VectorXd get_point(int i) const;
    MatrixXd collect() const;

    private:
    VectorXd start_point;
    VectorXd end_point;
    int grid_points_on_each_axis;
    int d;
};

#endif
