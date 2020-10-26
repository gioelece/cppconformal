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

## References

Zeni G, Fontana M, Vantini S. _Conformal Prediction: a Unified Review of Theory and New Challenges._ arXiv:200507972 [cs, econ, stat]. Published online May 16, 2020. Accessed October 26, 2020. http://arxiv.org/abs/2005.07972

Vovk V, Gammerman A, Shafer G. _Algorithmic Learning in a Random World._ Springer; 2005.

Nouretdinov I, Gammerman J, Fontana M, Rehal D. _Multi-level conformal clustering: A distribution-free technique for clustering and anomaly detection._ Neurocomputing. 2020;397:279-291. doi:10.1016/j.neucom.2019.07.114
