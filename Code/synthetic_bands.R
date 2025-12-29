
# ---- helper: make a synthetic SpatRaster with fixed N, variable B ----
make_synthetic_raster <- function(nrows = 500, ncols = 500, nbands = 3, seed = 1) {
  set.seed(seed)
  r <- rast(nrows = nrows, ncols = ncols, nlyrs = nbands,
            xmin = 0, xmax = ncols, ymin = 0, ymax = nrows)

  # Fill each band with random values
  # terra::values() can take a matrix with Ncells rows and nbands cols
  vals <- matrix(runif(ncell(r) * nbands), nrow = ncell(r), ncol = nbands)
  values(r) <- vals
  r
}

# ---- helper: approximate memory footprint (bytes) of outputs ----
# Focus on what you return: membership/distances (K layers each) + centers
result_bytes <- function(res) {
  as.numeric(object.size(res$distances)) +
    as.numeric(object.size(res$memberships)) +
    as.numeric(object.size(res$centers))
}

# ---- main experiment ----
band_scaling_experiment <- function(B_seq = 2:20,
                                   N_target = 2.5e5,
                                   num_clusters = 3,
                                   m = 2,
                                   seed = 1,
                                   do_plot = FALSE) {

  # Choose fixed dimensions matching N_target
  # Here: 500 x 500 = 250,000
  nrows <- 500
  ncols <- 500
  if (nrows * ncols != N_target) {
    stop("This script assumes N_target = 250000 (500x500). Adjust nrows/ncols if needed.")
  }

  results <- data.frame(
    B = integer(),
    N = integer(),
    K = integer(),
    runtime_sec = double(),
    output_MB = double()
  )

  for (B in B_seq) {
    img <- make_synthetic_raster(nrows = nrows, ncols = ncols, nbands = B, seed = seed)

    gc()
    t0 <- proc.time()[["elapsed"]]
    res <- im.fuzzy(
      input_image = img,
      num_clusters = num_clusters,
      seed = seed,
      m = m,
      do_plot = do_plot
    )
    runtime <- proc.time()[["elapsed"]] - t0

    out_mb <- result_bytes(res) / (1024^2)

    results <- rbind(
      results,
      data.frame(
        B = B,
        N = ncell(img),
        K = num_clusters,
        runtime_sec = runtime,
        output_MB = out_mb
      )
    )

    rm(img, res)
    gc()
  }

  results
}

# ---- run it ----
results_B <- band_scaling_experiment(B_seq = 2:20, N_target = 2.5e5, num_clusters = 4, m = 2, seed = 1)

# ---- corroborate "approximately linear runtime in B" ----
fit <- lm(runtime_sec ~ B, data = results_B)
print(summary(fit))

# ---- plot runtime vs bands ----
p_runtime_bands <- ggplot(results_B, aes(B, runtime_sec)) +
  geom_line(linewidth = 1, col = "green") +
  geom_point(size = 2.5, col = "green") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, col = "cyan1") +
  labs(
    title = "Runtime vs spectral dimensionality",
    subtitle = "Fixed image size N = 2.5 × 10^5 pixels\nFixed number of clusters K = 4",
    x = "Number of bands (B)",
    y = "Runtime (seconds)"
  ) +
  theme_minimal(base_size = 14)

print(p_runtime_bands)

