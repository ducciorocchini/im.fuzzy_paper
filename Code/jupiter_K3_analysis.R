# ============================================================
# JUPITER EMPIRICAL ANALYSIS - K = 3
#
# Outputs:
#   1. Fuzzy membership rasters
#   2. Membership map figure
#   3. Membership-weighted RGB brightness densities
#   4. Brightness density figure
#
# K = 3 was retained after evaluating candidate solutions
# from K = 2 to K = 10 using Silhouette, normalized FPC,
# and normalized Partition Entropy.
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(imageRy)
library(ggplot2)
library(ggridges)


# ------------------------------------------------------------
# 2. Import Jupiter image directly from GitHub
# ------------------------------------------------------------

jupiter_url <- paste0(
  "https://raw.githubusercontent.com/",
  "ducciorocchini/im.fuzzy_paper/main/Data/",
  "STScI-01EVSV9JYWYARC4A8336KMKNTQ.jpeg"
)

mg <- terra::rast(jupiter_url)

# Reproduce orientation used in the manuscript
mg <- terra::flip(mg)

print(mg)
terra::nlyr(mg)


# ------------------------------------------------------------
# 3. Optional display of the RGB image
# ------------------------------------------------------------

terra::plotRGB(
  mg,
  r = 1,
  g = 2,
  b = 3,
  stretch = "lin"
)


# ============================================================
# 4. Run im.fuzzy() with K = 3
# ============================================================

K <- 3
m <- 2
seed <- 42

res <- imageRy::im.fuzzy(
  mg,
  num_clusters = K,
  m = m,
  seed = seed,
  do_plot = FALSE
)


# ------------------------------------------------------------
# 5. Extract outputs
# ------------------------------------------------------------

mem <- res$memberships
dist <- res$distances
centers <- res$centers

print(mem)
print(dist)
print(centers)


# ============================================================
# 6. MEMBERSHIP MAPS
# ============================================================

names(mem) <- paste0(
  "Cluster ",
  seq_len(K)
)


# ------------------------------------------------------------
# Viridis-type continuous palette
# ------------------------------------------------------------

membership_palette <- hcl.colors(
  100,
  palette = "viridis"
)


# ------------------------------------------------------------
# Display membership maps
# ------------------------------------------------------------

terra::plot(
  mem,
  col = membership_palette,
  axes = FALSE,
  main = names(mem)
)


# ------------------------------------------------------------
# Save membership figure
# ------------------------------------------------------------

png(
  filename = "jupiter_memberships_K3.png",
  width = 2400,
  height = 1600,
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
# 7. RGB BRIGHTNESS
# ============================================================

# Mean RGB brightness for every pixel

brightness <- (
  mg[[1]] +
  mg[[2]] +
  mg[[3]]
) / 3

names(brightness) <- "brightness"


# ============================================================
# 8. RANDOM SAMPLE FOR FUZZY BRIGHTNESS DENSITIES
# ============================================================

set.seed(seed)

sample_size <- min(
  200000,
  terra::ncell(mg)
)

samp <- terra::spatSample(
  c(brightness, mem),
  size = sample_size,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

print(head(samp))


# ------------------------------------------------------------
# Separate brightness and membership matrix
# ------------------------------------------------------------

x <- samp$brightness

M <- as.matrix(
  samp[, -1, drop = FALSE]
)

colnames(M) <- paste0(
  "Cluster ",
  seq_len(K)
)


# ============================================================
# 9. Convert to long format
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
# 10. MEMBERSHIP-WEIGHTED BRIGHTNESS DENSITIES
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
    subtitle = "Jupiter image - K = 3 clusters",
    x = "RGB brightness",
    y = "Weighted density",
    fill = "Cluster"
  )

print(p_density)


# ------------------------------------------------------------
# Save brightness density figure
# ------------------------------------------------------------

ggsave(
  filename = "jupiter_fuzzy_brightness_K3.png",
  plot = p_density,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 11. OPTIONAL RIDGELINE VERSION
#
# This can be useful if overlapping density curves are
# difficult to distinguish.
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
    subtitle = "Jupiter image - K = 3 clusters",
    x = "RGB brightness",
    y = "Cluster"
  )

print(p_ridge)


ggsave(
  filename = "jupiter_fuzzy_brightness_ridges_K3.png",
  plot = p_ridge,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 12. SAVE NUMERICAL OUTPUTS
# ============================================================

write.csv(
  centers,
  "jupiter_centroids_K3.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# Save rasters if desired
# ------------------------------------------------------------

terra::writeRaster(
  mem,
  filename = "jupiter_memberships_K3.tif",
  overwrite = TRUE
)

terra::writeRaster(
  dist,
  filename = "jupiter_distances_K3.tif",
  overwrite = TRUE
)


# ============================================================
# 13. CHECK MEMBERSHIP SUM
#
# Memberships should sum to approximately 1 for each
# valid pixel.
# ============================================================

membership_sum <- sum(mem)

print(
  terra::global(
    membership_sum,
    c("min", "mean", "max"),
    na.rm = TRUE
  )
)

# output:
# jupiter_memberships_K3.png
# jupiter_fuzzy_brightness_K3.png
# jupiter_fuzzy_brightness_ridges_K3.png

