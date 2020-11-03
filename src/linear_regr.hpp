#ifndef __LINEAR_REGR_HPP
#define __LINEAR_REGR_HPP
#include <RcppEigen.h>
using namespace Eigen;

/*! Base class for linear regression models.
    Provides the default copy constructor.
*/
class LinearRegressionBase {
    public:
    /*! Use a fitted linear regression model to make a prediction
        \param X0 matrix of independent variables
        \returns prediction of Y corresponding to X0
    */
    template<typename Derived>
    MatrixXd predict(const MatrixBase<Derived> & X0) {
        if (!is_fitted) {
            Rcpp::stop("Linear model has not been fitted yet");
        }
        return X0 * beta;
    }

    protected:
    void set_beta(MatrixXd new_beta) {
        beta = new_beta;
        is_fitted = true;
    }

    private:
    MatrixXd beta;
    bool is_fitted = false;
};

/*! Class holding a linear regression model.
*/
class LinearRegression : public LinearRegressionBase {
    public:
    /*! Fit the linear regression model
        \param X matrix of independent variables
        \param Y matrix of covariates
    */
    template<typename Derived1, typename Derived2>
    void fit(const MatrixBase<Derived1> & X, const MatrixBase<Derived2> & y) {
        set_beta((X.transpose() * X).ldlt().solve(X.transpose() * y));
    }
};

/*! Class holding a ridge regression model.
*/
class RidgeRegression : public LinearRegressionBase {
    public:
    /*! Constructs a ridge regression model instance.
        \params l lambda (if lambda = 0, equivalent to linear regression)
    */
    RidgeRegression(double l) : lambda(l) {};

    /*! Fit the ridge regression model
        \param X matrix of independent variables
        \param Y matrix of covariates
    */
    template<typename Derived1, typename Derived2>
    void fit(const MatrixBase<Derived1> & X, const MatrixBase<Derived2> & Y) {
        int p = X.cols();
        auto eye = MatrixXd::Identity(p, p);
        set_beta((X.transpose() * X + lambda * eye).ldlt().solve(X.transpose() * Y));
    }

    private:
    double lambda;
};

#endif
