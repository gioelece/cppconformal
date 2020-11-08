library(ggplot2)
library(devtools)

# This loads the package in the current folder, without installing it
# (useful for development).
devtools::load_all()

Xhat = t(c(5, 1))
X = cbind(
    rnorm(200, sd=10),
    rnorm(200, sd=10)
)
sd = 0.5
y = cbind(
    X[, 1] + rnorm(200, sd=sd),
    2 * X[, 1] + rnorm(200, sd=sd)
)

res = run_linear_conformal_single_grid(as.matrix(X), y, as.matrix(Xhat))
# Or else, for example
# res = run_linear_conformal_multi_grid(
#     X, y, as.matrix(Xhat),
#     c(0.8, 0.9, 0.95), c(500, 500, 500, 500), 1.25
# )

grid = res$y_grid
p_values = res$p_values[1, ]

df = data.frame(
    y1 = grid[, 1], y2 = grid[, 2], p_values = p_values
)

df_large_p = df[df$p_values > 0.05, ]

ggplot(data = df, aes(x = y1, y = y2, z = p_values)) +
    geom_tile(aes(fill = p_values)) +
    xlim(min(df_large_p$y1) - sd, max(df_large_p$y1) + sd) +
    ylim(min(df_large_p$y2) - sd, max(df_large_p$y2) + sd) +
    stat_contour(breaks = c(0.05, 0.1)) +
    scale_fill_continuous()
ggsave("rnorm.example.png")
