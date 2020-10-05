#ifndef __LINEAR_REGR_HPP
#define __LINEAR_REGR_HPP
#include <RcppEigen.h>
using namespace Eigen;

class LinearRegression {
    public:
    void fit(MatrixXd const & X, MatrixXd const & y);
    MatrixXd predict(MatrixXd const & X0);

    private:
    MatrixXd beta;
}; 
#endif