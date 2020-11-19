/*! @file */
#ifndef __ALGORITHMS__BASE_HPP
#define __ALGORITHMS__BASE_HPP
#include <iostream>
#include <RcppEigen.h>
#include "../grid.hpp"

/*! Abstract class for a conformal algorithm
    \param Model class model base
*/
template<class Model>
class AlgorithmBase {
    public:
    /*! Run a conformal regression algorithm.

        __Remark__: multi_grid accepts only a single `Xhat`.

        \param model model to use as a base for conformal regression (will be copied at each run)
        \param X matrix of the independent variables
        \param Y matrix of the covariates
        \param Xhat a single point containing the values for the independent variables
        \return An Rcpp list
    */
    virtual List run(
        const Model & model,
        const MatrixXd & X, const MatrixXd & Y, const MatrixXd & Xhat
    ) = 0;
};

#endif
