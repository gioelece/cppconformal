#include "linear_regr.hpp"
using namespace Eigen;

void LinearRegression::fit(MatrixXd const & X, VectorXd const & y) {
    beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
}

VectorXd LinearRegression::predict(MatrixXd const & X0) {
    return X0 * beta;
};