#ifndef __GRID_HPP
#define __GRID_HPP
#include <vector>
#include <Eigen/Core>

using namespace Eigen;

class Grid {
    public:
    Grid(const VectorXd & s, const VectorXd & e, int g) :
        start_point(s), end_point(e), grid_side(g),
        d(s.size())
    {
        compute_step_increment();
    };

    const VectorXd & get_start_point() const {
        return start_point;
    };
    const VectorXd & get_end_point() const {
        return end_point;
    };
    const VectorXd & get_step_increment() const {
        return step_increment;
    };

    int get_size() const;
    VectorXd get_point(int point_idx) const;
    MatrixXd collect() const;

    private:
    void compute_step_increment();
    VectorXd start_point;
    VectorXd end_point;
    VectorXd step_increment;
    int grid_side;
    int d;
};

#endif
