#include "linear_regr.hpp"
using namespace Eigen;

void LinearRegression::fit(const MatrixXd & X, const MatrixXd & y) {
    beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
}

MatrixXd LinearRegression::predict(const MatrixXd & X0) {
    return X0 * beta;
};

void RidgeRegression::fit(const MatrixXd & X, const MatrixXd & y) {
    int p = X.cols();
    auto eye = MatrixXd::Identity(p, p);
    beta = (X.transpose() * X + lambda * eye).ldlt().solve(X.transpose() * y);
}

MatrixXd RidgeRegression::predict(const MatrixXd & X0) {
    return X0 * beta;
};
