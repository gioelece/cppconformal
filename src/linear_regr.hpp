#ifndef __LINEAR_REGR_HPP
#define __LINEAR_REGR_HPP
#include <RcppEigen.h>
using namespace Eigen;

class LinearRegression {
    public:
    template<typename Derived1, typename Derived2>
    void fit(const MatrixBase<Derived1> & X, const MatrixBase<Derived2> & y) {
        beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    }

    template<typename Derived>
    MatrixXd predict(const MatrixBase<Derived> & X0) {
        return X0 * beta;
    }

    private:
    MatrixXd beta;
};

class RidgeRegression {
    public:
    RidgeRegression(double l) : lambda(l) {};

    template<typename Derived1, typename Derived2>
    void fit(const MatrixBase<Derived1> & X, const MatrixBase<Derived2> & y) {
        int p = X.cols();
        auto eye = MatrixXd::Identity(p, p);
        beta = (X.transpose() * X + lambda * eye).ldlt().solve(X.transpose() * y);
    }

    template<typename Derived>
    MatrixXd predict(const MatrixBase<Derived> & X0) {
        return X0 * beta;
    }

    private:
    double lambda;
    MatrixXd beta;
};

#endif
