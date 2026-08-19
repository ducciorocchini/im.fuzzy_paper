# ============================================================
# FINAL JUPITER ANALYSIS - MASKED IMAGE, K = 2
# ============================================================

library(terra)
library(imageRy)
library(ggplot2)
library(ggridges)

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------

K <- 2
m <- 2
seed <- 42

# ------------------------------------------------------------
# 2. Run im.fuzzy() on masked Jupiter image
# ------------------------------------------------------------

res <- imageRy::im.fuzzy(
  img_masked,
  num_clusters = K,
  m = m,
  seed = seed,
  do_plot = FALSE
)

mem <- res$memberships
dist <- res$distances
centers <- res$centers

names(mem) <- paste0("Cluster ", seq_len(K))

print(mem)
print(centers)

# ============================================================
# 3. FINAL MEMBERSHIP MAPS
# ============================================================

membership_palette <- hcl.colors(
  100,
  palette = "viridis"
)

terra::plot(
  mem,
  col = membership_palette,
  axes = FALSE,
  main = names(mem)
)

png(
  filename = "jupiter_memberships_masked_K2.png",
  width = 2400,
  height = 1300,
  res = 300
)

terra::plot(
  mem,
  col = membership_palette,
  axes = FALSE,
  main = names(mem)
)

dev.off()

# ============================================================
# 4. RGB brightness on masked image
# ============================================================

brightness <- (
  img_masked[[1]] +
  img_masked[[2]] +
  img_masked[[3]]
) / 3

names(brightness) <- "brightness"

# ============================================================
# 5. Sample pixels for brightness distributions
# ============================================================

set.seed(seed)

sample_size <- min(
  200000,
  terra::ncell(img_masked)
)

samp <- terra::spatSample(
  c(brightness, mem),
  size = sample_size,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

x <- samp$brightness

M <- as.matrix(
  samp[, -1, drop = FALSE]
)

colnames(M) <- paste0(
  "Cluster ",
  seq_len(K)
)

# ============================================================
# 6. Long format
# ============================================================

long <- data.frame(
  brightness = rep(
    x,
    times = K
  ),
  cluster = factor(
    rep(
      colnames(M),
      each = length(x)
    ),
    levels = colnames(M)
  ),
  weight = as.vector(M)
)

# ============================================================
# 7. Ridgeline plot
# ============================================================

p_ridge <- ggplot(
  long,
  aes(
    x = brightness,
    y = cluster,
    weight = weight,
    fill = cluster
  )
) +
  ggridges::geom_density_ridges(
    alpha = 0.7,
    scale = 1.1,
    rel_min_height = 0.01
  ) +
  scale_fill_viridis_d(
    option = "viridis",
    end = 0.9
  ) +
  theme_minimal(
    base_size = 14
  ) +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "Membership-weighted brightness distributions",
    subtitle = "Masked Jupiter image - K = 2 clusters",
    x = "RGB brightness",
    y = "Cluster"
  )

print(p_ridge)

ggsave(
  filename = "jupiter_fuzzy_brightness_ridges_masked_K2.png",
  plot = p_ridge,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# 8. Optional overlapping density version
# ============================================================

p_density <- ggplot(
  long,
  aes(
    x = brightness,
    weight = weight,
    fill = cluster
  )
) +
  geom_density(
    alpha = 0.45,
    adjust = 1.2
  ) +
  scale_fill_viridis_d(
    option = "viridis",
    end = 0.9
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Fuzzy brightness densities",
    subtitle = "Masked Jupiter image - K = 2 clusters",
    x = "RGB brightness",
    y = "Weighted density",
    fill = "Cluster"
  )

print(p_density)

ggsave(
  filename = "jupiter_fuzzy_brightness_masked_K2.png",
  plot = p_density,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# 9. Save centroids and rasters
# ============================================================

write.csv(
  centers,
  "jupiter_centroids_masked_K2.csv",
  row.names = FALSE
)

terra::writeRaster(
  mem,
  "jupiter_memberships_masked_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  dist,
  "jupiter_distances_masked_K2.tif",
  overwrite = TRUE
)

# ============================================================
# 10. Check membership sum
# ============================================================

membership_sum <- sum(mem)

print(
  terra::global(
    membership_sum,
    c("min", "mean", "max"),
    na.rm = TRUE
  )
)
