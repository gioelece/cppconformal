#ifndef __GRID_HPP
#define __GRID_HPP
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;
using Rcpp::List;

/*! Class holding a (hyper-)rectangular grid, that avoids storing in memory the coordinates of each point.
*/  
class Grid {
    public:
    /*! Construct a grid object.
        \param s start point (bottom-left)
        \param e end point (top-right)
        \param g number of points for each side of the grid
        \return The point coordinates 
    */  
    Grid(const VectorXd & s, const VectorXd & e, int g) :
        start_point(s), end_point(e), grid_side(g),
        d(s.size())
    {
        compute_step_increment();
    };

    /*! Get the starting point of the grid (bottom-left).
        \return The point coordinates 
    */  
    const VectorXd & get_start_point() const {
        return start_point;
    };
    
    /*! Get the end point of the grid (top-right).
        \return The point coordinates 
    */  
    const VectorXd & get_end_point() const {
        return end_point;
    };

    /*! Get the length of a grid step.
        \return A vector with grid steps length (one for each direction)
    */  
    const VectorXd & get_step_increment() const {
        return step_increment;
    };

    /*! Get number of points for each side of the grid.
    */  
    int get_grid_side() const {
        return grid_side;
    };

    /*! Get the size of the grid.
        \return The number of grid points for this grid
    */  
    int get_size() const;

    /*! Get the i-th point of the grid (0: bottom left, get_size() - 1: top right), where i = point_idx;
        This is not stored, but calcolated on the fly.
        \return The point coordinates 
    */  
    VectorXd get_point(int point_idx) const;

    /*! Compute the coordinates for each point of the grid and collect them in a matrix.
        \return The point coordinates 
    */  
    MatrixXd collect() const;

    /*! Get the parameters of the grid as an Rcpp list
        \return A Rcpp list
    */
    List get_parameters_as_list() const;

    private:
    void compute_step_increment();
    VectorXd start_point;
    VectorXd end_point;
    VectorXd step_increment;
    int grid_side;
    int d;
};

#endif
