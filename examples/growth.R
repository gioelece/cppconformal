library(fda)
library(ggplot2)
library(devtools)

# This loads the package in the current folder, without installing it
# (useful for development).
devtools::load_all()

p = 5
n_boys = dim(growth$hgtm)[2]
n_girls = dim(growth$hgtf)[2]

# Columns of X, X_boys, X_girls:
# 1: is_girl, 2: is_boy
X_boys = matrix(rep(c(0, 1), n_boys), ncol=2, byrow=TRUE)
X_girls = matrix(rep(c(1, 0), n_girls), ncol=2, byrow=TRUE)

age_selection = c(20, 21, 22)
Y_boys = t(growth$hgtm[age_selection, ])
Y_girls = t(growth$hgtf[age_selection, ])

X = rbind(X_boys, X_girls)
Y = rbind(Y_boys, Y_girls)

X0 = rbind(c(0, 1)) # Boys

res = run_linear_conformal_multi_grid(
    X, Y, X0,
    c(0.7, 0.8, 0.9, 0.95), c(20, 10, 10, 10, 100), 1.00
)

grid = res$y_grid
p_values = res$p_values[1, ]
df = data.frame(
    y1 = grid[, 1], y2 = grid[, 2], p_values = p_values
)

df_large_p = df[df$p_values > 0.8, ]

ggplot(data = df_large_p, aes(x = y1, y = y2, z = p_values)) +
    geom_tile(aes(fill = p_values)) +
    # xlim(min(df_large_p$y1), max(df_large_p$y1)) +
    # ylim(min(df_large_p$y2), max(df_large_p$y2)) +
    stat_contour() +
    scale_fill_continuous()
ggsave("growth.example.png")
