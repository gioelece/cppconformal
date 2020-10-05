# cppconformal

An implementation of conformal regression with Rcpp.

## How to install
```r
library('devtools')
devtools::install_github('gioelece/cppconformal')
```

## How to use

One can call `run_linear_conformal(X, Y, X0, grid_size=150, grid_param=1.25)`:
where `X` ($n \times p$ matrix) contains the covariates, `Y` ($n \times d$) is the matrix of the corresponding observed values, and one wants to construct a confidence interval for the response to the covariates `X0` ($n_0 \times p$).

To sample the response space, a uniform grid is created. The limits of the grid for the $i$-th axis are `-limit_i` to `+limit_i` where `limit_i = grid_param * max(abs(y_i))`, with `grid_size` points.

**Remark**: the intercept coefficient is not included in the prediction. To have a "typical" linear regression, one needs to add to `X` a column of ones.