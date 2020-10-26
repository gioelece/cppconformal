# cppconformal

An implementation of conformal regression with Rcpp.

## How to install
```r
library('devtools')
devtools::install_github('gioelece/cppconformal')
```

## How to use

The library exports in R the following functions:
```R
run_linear_conformal_single_grid(X, y, X0, grid_side, grid_param)
run_ridge_conformal_single_grid(X, y, X0, lambda, grid_side, grid_param)
run_linear_conformal_multi_grid(X, y, X0, grid_levels, grid_sides, initial_grid_param)
run_ridge_conformal_multi_grid(X, y, X0, lambda, grid_levels, grid_sides, initial_grid_param)
```

For example, one can call `run_linear_conformal(X, Y, X0, grid_side, grid_param)`:
`X` ($n \times p$ matrix) contains the covariates, `Y` ($n \times d$) is the matrix of the corresponding observed values, and one wants to construct a confidence interval for the response to the covariates `X0` ($n_0 \times p$).

To sample the response space, for the `single_grid` function family, a uniform grid is created. The limits of the grid for the $i$-th axis are `-limit_i` to `+limit_i` where `limit_i = grid_param * max(abs(y_i))`, with `grid_side` points for each dimension.

Instead, when using a `*_multi_grid` function, an initial "coarse" grid is created as before, with parameters `initial_grid_param` and `grid_sides[0]`. Then a subgrid of size `grid_sides[1]` is created to contain all the points (from the previous grid) where the value of $p$ is greater or equal than `grid_levels[0]`, and so on, for all the elements of `grid_levels`. Note that, in order to use these functions, one needs to have a single `X0`, i.e. $n_0 = 1$.

Let $g = grid_side ^ d$ be the total number of grid points. The functions return a R list with `grid` ($g \times d$), containing the sampled points, and `p_values` ($n_0 \times g$), containing the corresponding p-values for each `X0`. For `*_multi_grid` functions, only the values referring to the last grid are returned.

**Remark**: the intercept coefficient is not included in the prediction. To have a "typical" linear regression, one needs to add to `X` a column of ones.

## References

Zeni G, Fontana M, Vantini S. _Conformal Prediction: a Unified Review of Theory and New Challenges._ arXiv:200507972 [cs, econ, stat]. Published online May 16, 2020. Accessed October 26, 2020. http://arxiv.org/abs/2005.07972

Vovk V, Gammerman A, Shafer G. _Algorithmic Learning in a Random World._ Springer; 2005.

Nouretdinov I, Gammerman J, Fontana M, Rehal D. _Multi-level conformal clustering: A distribution-free technique for clustering and anomaly detection._ Neurocomputing. 2020;397:279-291. doi:10.1016/j.neucom.2019.07.114
