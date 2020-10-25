library(ggplot2)
library(devtools)

# This loads the package in the current folder, without installing it
# (useful for development).
devtools::load_all()

X0 = 5
X = rnorm(200, sd=10)
y = cbind(
    X + rnorm(200, sd=1),
    2 * X + rnorm(200, sd=1)
)

# res = run_linear_conformal_single_grid(as.matrix(X), y, as.matrix(X0))
res = run_linear_conformal_multi_grid(
    as.matrix(X), y, as.matrix(X0),
    c(0.8, 0.9, 0.95), c(500, 500, 500, 500), 1.25
)

grid = res$y_grid
p_values = res$p_values[1, ]

df = data.frame(
    y1 = grid[, 1], y2 = grid[, 2], p_values = p_values
)

df_large_p = df[df$p_values > 0.5, ]

ggplot(data = df, aes(x = y1, y = y2, z = p_values)) +
    geom_tile(aes(fill = p_values)) +
    xlim(min(df_large_p$y1), max(df_large_p$y1)) +
    ylim(min(df_large_p$y2), max(df_large_p$y2)) +
    stat_contour() +
    scale_fill_continuous()
ggsave("example.png")
