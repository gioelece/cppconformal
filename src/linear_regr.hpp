#ifndef __LINEAR_REGR_HPP
#define __LINEAR_REGR_HPP
#include <RcppEigen.h>
using namespace Eigen;

class LinearRegression {
    public:
    void fit(const MatrixXd & X, const MatrixXd & y);
    MatrixXd predict(const MatrixXd & X0);

    private:
    MatrixXd beta;
};

class RidgeRegression {
    public:
    RidgeRegression(double l) : lambda(l) {};
    void fit(const MatrixXd & X, const MatrixXd & y);
    MatrixXd predict(const MatrixXd & X0);

    private:
    double lambda;
    MatrixXd beta;
};

#endif
