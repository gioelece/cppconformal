library(fda)
library(ggplot2)
library(reshape)
library(devtools)

# This loads the package in the current folder, without installing it
# (useful for development).
# devtools::load_all()
library(cppconformal)

n_boys = dim(growth$hgtm)[2]
n_girls = dim(growth$hgtf)[2]
ages = rownames(growth$hgtf)

# Columns of X, X_boys, X_girls:
# 1: is_girl, 2: is_boy
X_boys = matrix(rep(c(0, 1), n_boys), ncol=2, byrow=TRUE)
X_girls = matrix(rep(c(1, 0), n_girls), ncol=2, byrow=TRUE)

age_selection = c(1, 6, 11, 16, 21, 26, 31)
Y_boys = t(growth$hgtm[age_selection, ])
Y_girls = t(growth$hgtf[age_selection, ])

X = rbind(X_boys, X_girls)
Y = as.data.frame(rbind(Y_boys, Y_girls))
Yscaled = scale(Y)

Xhat = rbind(c(0, 1)) # Boys

ptm = proc.time()
res = run_linear_conformal_multi_grid(
    X, Yscaled, Xhat,
    c(0.01, 0.05), c(10, 10, 10), 3,
    print_progress = TRUE
)
elapsed = proc.time() - ptm
elapsed

data_filter = order(res$p_values[1, ], decreasing=TRUE)[1:1000]
yscaled_pts = as.data.frame(res$y_grid[data_filter, ])
y_pts = yscaled_pts * attr(Yscaled, "scaled:scale") + attr(Yscaled, "scaled:center")
data = as.data.frame(cbind(y_pts, res$p_values[1, data_filter]))
colnames(data) <- c(ages[age_selection], "p_value")

plot_data = data #[sample(nrow(data), size=2000, prob=data$p_value), ]
sel_ages = ages[age_selection]
new_grid = plot_data[, 1:length(age_selection)]
new_p_values = plot_data[, length(age_selection) + 1]


# Plots
png('growth_dataset.example.png', width=700, height=700)
matplot(ages, growth$hgtm, type='l', xlab="Age [years]", ylab="Height [cm]",
    col=rgb(red=0, green=0, blue=1, alpha=0.3),
    cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
matplot(ages, growth$hgtf, type='l', add=TRUE,
    col=rgb(red=1, green=0, blue=0, alpha=0.3))
dev.off()

png('growth.example.png', width=700, height=700)
matplot(sel_ages, t(new_grid), type='l',
    col=rgb(red=0, green=0, blue=1, alpha=new_p_values * 0.1),
    cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
dev.off()

png('growth_grid.example.png', width=700, height=700)
grid_bounds = t(sapply(res$y_grid_parameters, function(g) {
    scaled_start = g$start_point * attr(Yscaled, "scaled:scale") + attr(Yscaled, "scaled:center")
    scaled_end = g$end_point * attr(Yscaled, "scaled:scale") + attr(Yscaled, "scaled:center")
    c(scaled_start[1:2], scaled_end[1:2])
}))

plot(cbind(grid_bounds[, 1], grid_bounds[, 3]),
     cbind(grid_bounds[, 2], grid_bounds[, 4]), type="n",
     xlab="Y_1", ylab="Y_2",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
rect(grid_bounds[, 1], grid_bounds[, 2], grid_bounds[, 3], grid_bounds[, 4])
dev.off()

# get_prediction_for_age = function(i) {
#     res = run_linear_conformal_multi_grid(
#         X, as.matrix(Y[, i]), Xhat,
#         c(0.001), c(10, 2000), 1.00,
#     )
#     range(res$y_grid[res$p_values > 0.05 / 7, ])
# }

# bonferroni_pred = sapply(seq_len(length(age_selection)), get_prediction_for_age)
# lines(sel_ages, bonferroni_pred[1, ])
# lines(sel_ages, bonferroni_pred[2, ])
