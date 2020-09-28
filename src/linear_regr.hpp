#ifndef __LINEAR_REGR_HPP
#define __LINEAR_REGR_HPP
#include <RcppEigen.h>
using namespace Eigen;

class LinearRegression {
    public:
    void fit(MatrixXd const & X, VectorXd const & y);
    VectorXd predict(MatrixXd const & X0);

    private:
    VectorXd beta;
}; 
#endif