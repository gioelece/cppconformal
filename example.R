# This loads the package in the current folder, without installing it
# (useful for development).
require("devtools")
devtools::load_all()

X = rnorm(200, sd=10)
X0 = rnorm(1, sd=10)
y = X + rnorm(200, sd=1)

X0
run_linear_conformal(as.matrix(X), y, as.matrix(X0))