# --- Packages ---------------------------------------------------------------
library(terra)
library(ggplot2)
library(imageRy)   # contains im.fuzzy()

# --- Synthetic raster (fixed N, fixed B) ------------------------------------
make_synthetic_raster <- function(n = 500, B = 4, seed = 1) {
  set.seed(seed)
  r <- rast(nrows = n, ncols = n, nlyrs = B)
  values(r) <- matrix(runif(ncell(r) * B), nrow = ncell(r), ncol = B)
  names(r) <- paste0("band", 1:B)
  r
}

# --- Experiment -------------------------------------------------------------
cluster_scaling_experiment <- function(
  K_seq = 2:12,
  n = 500,          # 500 x 500 = 250k pixels
  B = 4,
  m = 2,
  seed = 1
) {

  img <- make_synthetic_raster(n = n, B = B, seed = seed)
  N <- ncell(img)

  results <- data.frame(
    K = integer(),
    N = integer(),
    B = integer(),
    runtime_sec = double()
  )

  for (K in K_seq) {

    cat("Running K =", K, "\n")
    gc()

    t0 <- system.time({
      out <- im.fuzzy(
        img,
        num_clusters = K,
        m = m,
        do_plot = FALSE
      )
    })[["elapsed"]]

    results <- rbind(
      results,
      data.frame(
        K = K,
        N = N,
        B = B,
        runtime_sec = t0
      )
    )

    rm(out)
    gc()
  }

  results
}

# --- Run experiment ---------------------------------------------------------
results_K <- cluster_scaling_experiment(
  K_seq = 2:12,
  n = 500,
  B = 4,
  m = 2,
  seed = 1
)

print(results_K)

# --- Plot -------------------------------------------------------------------
p_runtime_clusters <- ggplot(results_K, aes(x = K, y = runtime_sec)) +
  geom_line(linewidth = 1, color = "green") +
  geom_point(size = 3, color = "green") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, color = "cyan1") +
  scale_x_continuous(
    breaks = seq(2, max(results_K$K), by = 2),     # ⬅ grid lines at 2, 4, 6, …
    limits = c(min(results_K$K), max(results_K$K))
  ) +
  labs(
    title = "Runtime vs number of clusters",
    subtitle = "Fixed image size N = 2.5 × 10⁵ pixels\nFixed number of bands B = 4",
    x = "Number of clusters (K)",
    y = "Runtime (seconds)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor.x = element_blank()           # ⬅ no odd/fractional grids
  )

print(p_runtime_clusters)

