
# --- Packages ---------------------------------------------------------------
library(terra)
library(imageRy)   # contains im.fuzzy()

# --- Generate synthetic multiband raster ----------------------------
make_synthetic_raster <- function(n, B = 4, K = 3, overlap = 0.3, seed = 100) {
  set.seed(seed)

  # grid
  r <- rast(nrows = n, ncols = n, nlyrs = B)

  # class means
  means <- matrix(
    runif(K * B, min = 0, max = 100),
    nrow = K, ncol = B
  )

  # class covariance
  Sigma <- diag(B) * (1 + overlap)

  # assign latent class per pixel
  cls <- sample(1:K, ncell(r), replace = TRUE)

  # generate pixel values
  vals <- matrix(NA, nrow = ncell(r), ncol = B)
  for (k in 1:K) {
    idx <- which(cls == k)
    if (length(idx) > 0) {
      vals[idx, ] <- MASS::mvrnorm(length(idx), means[k, ], Sigma)
    }
  }

  values(r) <- vals
  names(r) <- paste0("band", 1:B)
  r
}

# --- Experimental settings --------------------------------------------------
sizes <- c(100, 250, 500, 750, 1000)  # image side length
B <- 4
K <- 3
m <- 2

results <- data.frame(
  n = sizes,
  N = sizes^2,
  runtime_sec = NA
)

# --- Run experiments --------------------------------------------------------
for (i in seq_along(sizes)) {

  cat("Running size:", sizes[i], "x", sizes[i], "\n")

  img <- make_synthetic_raster(
    n = sizes[i],
    B = B,
    K = K,
    seed = 100 + i
  )

  t0 <- Sys.time()

  out <- im.fuzzy(
    img,
    num_clusters = K,
    m = 2
  )

  t1 <- Sys.time()

  results$runtime_sec[i] <- as.numeric(difftime(t1, t0, units = "secs"))

  rm(img, out)
  gc()
}

print(results)

# --- Simple scaling plot ------------------------------------------
p_runtime_size <- ggplot(results, aes(x = N, y = runtime_sec)) +
  geom_line(linewidth = 1, col = "green") +
  geom_point(size = 3, col = "green") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, col = "cyan1") +
  labs(
    title = "Runtime vs image size",
    subtitle = "Fixed number of bands B = 4\nFixed number of clusters K = 4",
    x = "Number of pixels (N)",
    y = "Runtime (seconds)"
  ) +
  theme_minimal(base_size = 14)

print(p_runtime_size)


