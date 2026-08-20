# ============================================================
# MONSUMMANO - FUZZY MEMBERSHIP COMPARISON
#
# Panels:
# (a) Original RGB image
# (b) k-means crisp partition
# (c) fuzzy C-means membership - Cluster 1
# (d) fuzzy C-means membership - Cluster 2
# (e) im.fuzzy() membership - Cluster 1
# (f) im.fuzzy() membership - Cluster 2
#
# K = 2
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(e1071)
library(imageRy)
library(clue)
library(viridisLite)


# ------------------------------------------------------------
# 2. Load Monsummano image
# ------------------------------------------------------------

monsummano_url <- paste0(
  "https://raw.githubusercontent.com/",
  "ducciorocchini/im.fuzzy_paper/main/Data/",
  "monsummano.png"
)

img <- terra::rast(
  monsummano_url
)

img <- flip(img)

print(img)


# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------

K <- 2
m <- 2
seed <- 42


# ------------------------------------------------------------
# 4. Extract valid pixels
# ------------------------------------------------------------

X <- terra::as.matrix(
  img
)

valid_idx <- complete.cases(
  X
)

X_valid <- X[
  valid_idx,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Scale each band to 0-255
# ------------------------------------------------------------

scale_255 <- function(x) {

  xmin <- min(
    x,
    na.rm = TRUE
  )

  xmax <- max(
    x,
    na.rm = TRUE
  )

  if (xmax == xmin) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  255 * (x - xmin) / (xmax - xmin)
}


X_scaled <- apply(
  X_valid,
  2,
  scale_255
)

X_scaled <- as.matrix(
  X_scaled
)


# ============================================================
# 6. K-MEANS
# ============================================================

set.seed(seed)

km <- stats::kmeans(
  X_scaled,
  centers = K
)

labels_km <- km$cluster


# ============================================================
# 7. FUZZY C-MEANS
# ============================================================

set.seed(seed)

fcm <- e1071::cmeans(
  X_scaled,
  centers = K,
  m = m
)

U_fcm <- fcm$membership

labels_fcm <- apply(
  U_fcm,
  1,
  which.max
)


# ============================================================
# 8. im.fuzzy()
# ============================================================

imf <- imageRy::im.fuzzy(
  img,
  num_clusters = K,
  m = m,
  seed = seed,
  do_plot = FALSE
)

U_imf_all <- terra::values(
  imf$memberships,
  mat = TRUE
)

U_imf <- U_imf_all[
  valid_idx,
  ,
  drop = FALSE
]

labels_imf <- apply(
  U_imf,
  1,
  which.max
)


# ============================================================
# 9. MATCH CLUSTER IDENTITIES TO K-MEANS
# ============================================================

match_memberships <- function(
  reference_labels,
  target_labels,
  U_target,
  K
) {

  tab <- table(
    factor(
      reference_labels,
      levels = seq_len(K)
    ),
    factor(
      target_labels,
      levels = seq_len(K)
    )
  )

  assignment <- clue::solve_LSAP(
    tab,
    maximum = TRUE
  )

  # For reference cluster 1,2,... find the corresponding
  # target membership column
  matched_columns <- as.integer(
    assignment
  )

  U_target[
    ,
    matched_columns,
    drop = FALSE
  ]
}


U_fcm_matched <- match_memberships(
  reference_labels = labels_km,
  target_labels = labels_fcm,
  U_target = U_fcm,
  K = K
)


U_imf_matched <- match_memberships(
  reference_labels = labels_km,
  target_labels = labels_imf,
  U_target = U_imf,
  K = K
)


# ============================================================
# 10. BUILD K-MEANS PARTITION RASTER
# ============================================================

template <- img[[1]]

r_km <- template

vals_km <- rep(
  NA_integer_,
  terra::ncell(template)
)

vals_km[valid_idx] <- labels_km

terra::values(
  r_km
) <- vals_km


# ============================================================
# 11. HELPER FUNCTION FOR MEMBERSHIP RASTERS
# ============================================================

make_membership_raster <- function(
  membership_values
) {

  r <- template

  vals <- rep(
    NA_real_,
    terra::ncell(template)
  )

  vals[valid_idx] <- membership_values

  terra::values(r) <- vals

  r
}


# ============================================================
# 12. BUILD ALL FOUR FUZZY MEMBERSHIP RASTERS
# ============================================================

r_fcm_1 <- make_membership_raster(
  U_fcm_matched[, 1]
)

r_fcm_2 <- make_membership_raster(
  U_fcm_matched[, 2]
)

r_imf_1 <- make_membership_raster(
  U_imf_matched[, 1]
)

r_imf_2 <- make_membership_raster(
  U_imf_matched[, 2]
)


# ============================================================
# 13. CHECK MEMBERSHIP SUMS
# ============================================================

cat(
  "\nFCM membership sum range:\n"
)

print(
  range(
    rowSums(U_fcm_matched),
    na.rm = TRUE
  )
)

cat(
  "\nim.fuzzy membership sum range:\n"
)

print(
  range(
    rowSums(U_imf_matched),
    na.rm = TRUE
  )
)


# ============================================================
# 14. NUMERICAL MEMBERSHIP DIFFERENCE
# ============================================================

difference_matrix <- abs(
  U_fcm_matched -
  U_imf_matched
)

cat("\n")
cat("============================================\n")
cat("FUZZY MEMBERSHIP DIFFERENCE\n")
cat("============================================\n")

cat(
  "Mean absolute difference:",
  mean(
    difference_matrix,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median absolute difference:",
  median(
    difference_matrix,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum absolute difference:",
  max(
    difference_matrix,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 15. COLORS
# ============================================================

cluster_colors <- viridisLite::viridis(
  K,
  option = "D",
  end = 0.9
)

membership_colors <- viridisLite::viridis(
  100,
  option = "D"
)


# ============================================================
# 16. OUTPUT MAP FIGURES
# ============================================================

# ------------------------------------------------------------
# (a) RGB image
# ------------------------------------------------------------

imageRy::im.plotRGB(
  img,
  r = 1,
  g = 2,
  b = 3#,
  #title = "(a) UAV RGB image"
)

# ------------------------------------------------------------
# k-means partition
# ------------------------------------------------------------

terra::plot(
  r_km,
  col = cluster_colors,
  breaks = seq(
    0.5,
    K + 0.5,
    by = 1
  ),
  axes = T,
  legend = FALSE,
  main = "k-means partition"
)

# Get raster extent
e <- terra::ext(r_km)

# Put legend safely inside upper-right corner
legend(
  x = e$xmax - 0.01 * (e$xmax - e$xmin),
  y = e$ymax - 0.01 * (e$ymax - e$ymin),
  legend = c(
    "Cluster 1",
    "Cluster 2"
  ),
  fill = cluster_colors,
  border = "black",
  bg = "white",
  bty = "o",
  cex = 0.8,
  xjust = 1,
  yjust = 1
)

# ------------------------------------------------------------
# (a) FCM Cluster 1
# ------------------------------------------------------------

png(
  filename = "fuzzy_c_means.png",
  width = 3200,
  height = 4200,
  res = 300
)

par(
  mfrow = c(1, 2),
  mar = c(1, 1, 3, 3)
)

terra::plot(
  r_fcm_1,
  col = membership_colors,
  range = c(0, 1),
  axes = T,
  main = "(a) fuzzy C-means - Cluster 1"
)


# ------------------------------------------------------------
# (b) FCM Cluster 2
# ------------------------------------------------------------

terra::plot(
  r_fcm_2,
  col = membership_colors,
  range = c(0, 1),
  axes = T,
  main = "(b) fuzzy C-means - Cluster 2"
)

dev.off()

# ============================================================
# 17. SAVE MEMBERSHIP RASTERS
# ============================================================

terra::writeRaster(
  r_fcm_1,
  "monsummano_fcm_membership_cluster1_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_fcm_2,
  "monsummano_fcm_membership_cluster2_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_imf_1,
  "monsummano_imfuzzy_membership_cluster1_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_imf_2,
  "monsummano_imfuzzy_membership_cluster2_K2.tif",
  overwrite = TRUE
)

##########
##########

# ============================================================
# MONSUMMANO
# MEMBERSHIP-WEIGHTED RGB BRIGHTNESS DISTRIBUTIONS
# ============================================================

library(ggplot2)
library(viridisLite)


# ------------------------------------------------------------
# 1. RGB brightness
#
# Arithmetic mean of R, G and B for each valid pixel
# ------------------------------------------------------------

brightness_all <- rowMeans(
  X_valid[, 1:3, drop = FALSE]
)


# ------------------------------------------------------------
# 2. Build long data frame
#
# Use the matched im.fuzzy memberships so Cluster 1 and
# Cluster 2 correspond to those shown in the membership maps.
# ------------------------------------------------------------

brightness_long <- rbind(

  data.frame(
    brightness = brightness_all,
    Cluster = "Cluster 1",
    Membership = U_imf_matched[, 1]
  ),

  data.frame(
    brightness = brightness_all,
    Cluster = "Cluster 2",
    Membership = U_imf_matched[, 2]
  )
)


brightness_long$Cluster <- factor(
  brightness_long$Cluster,
  levels = c(
    "Cluster 1",
    "Cluster 2"
  )
)


# ------------------------------------------------------------
# 3. Same purple / lime palette used in the paper
# ------------------------------------------------------------

paper_colors <- viridisLite::viridis(
  2,
  option = "D",
  end = 0.9
)


# ============================================================
# 4. MEMBERSHIP-WEIGHTED DENSITY PLOT
# ============================================================

p_brightness <- ggplot(
  brightness_long,
  aes(
    x = brightness,
    weight = Membership,
    fill = Cluster
  )
) +

  geom_density(
    alpha = 0.55,
    adjust = 1.2,
    color = "black",
    linewidth = 0.8
  ) +

  scale_fill_manual(
    values = paper_colors
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Fuzzy brightness densities",
    subtitle = "Monsummano image - K = 2 clusters",
    x = "RGB brightness",
    y = "Weighted density",
    fill = "Cluster"
  ) +

  theme(
    legend.position = "right"
  )


# ============================================================
# 5. DISPLAY
# ============================================================

print(
  p_brightness
)


# ============================================================
# 6. SAVE
# ============================================================

ggsave(
  filename = "monsummano_fuzzy_brightness_K2.png",
  plot = p_brightness,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# im.fuzzy() UNCERTAINTY MAP
#
# uncertainty = 1 - maximum membership
#
# For K = 2:
#   0.0 = high certainty
#   0.5 = maximum uncertainty
# ============================================================

uncertainty_imf <- 1 - apply(
  U_imf_matched,
  1,
  max
)

# Build raster
r_imf_uncertainty <- template

uncertainty_vals <- rep(
  NA_real_,
  terra::ncell(template)
)

uncertainty_vals[valid_idx] <- uncertainty_imf

terra::values(
  r_imf_uncertainty
) <- uncertainty_vals


# ------------------------------------------------------------
# Check range
# ------------------------------------------------------------

cat("\n")
cat("im.fuzzy uncertainty range:\n")

print(
  range(
    uncertainty_imf,
    na.rm = TRUE
  )
)


#### Export paper maps

png(
  filename = "monsummano_fuzzy_memberships_K2.png",
  width = 3200,
  height = 4200,
  res = 300
)

par(
  mfrow = c(1, 3),
  mar = c(1, 1, 3, 3)
)


# ------------------------------------------------------------
# (e) im.fuzzy Cluster 1
# ------------------------------------------------------------

terra::plot(
  r_imf_1,
  col = membership_colors,
  range = c(0, 1),
  axes = T,
  main = "(a) im.fuzzy() - Cluster 1"
)


# ------------------------------------------------------------
# (f) im.fuzzy Cluster 2
# ------------------------------------------------------------

terra::plot(
  r_imf_2,
  col = membership_colors,
  range = c(0, 1),
  axes = T,
  main = "(b) im.fuzzy() - Cluster 2"
)

# plot
terra::plot(r_imf_uncertainty, col=mako(100), main = "(c) Membership uncertainty")

dev.off()
