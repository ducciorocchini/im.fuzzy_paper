# ============================================================
# MONSUMMANO - FUZZY PARTITION COMPARISON
#
# Panels:
# (a) Original RGB image
# (b) k-means crisp partition
# (c) fuzzy C-means membership - Cluster 1
# (d) im.fuzzy() membership - Cluster 1
# (e) absolute membership difference
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
#
# Same preprocessing used internally by im.fuzzy()
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

set.seed(
  seed
)

km <- stats::kmeans(
  X_scaled,
  centers = K
)

labels_km <- km$cluster


# ============================================================
# 7. FUZZY C-MEANS
# ============================================================

set.seed(
  seed
)

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


# ------------------------------------------------------------
# Extract membership matrix
# ------------------------------------------------------------

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
# 9. MATCH FCM AND im.fuzzy CLUSTERS TO K-MEANS
#
# Cluster IDs are arbitrary.
# We use k-means as the reference so that "Cluster 1"
# corresponds to the same radiometric group in all methods.
# ============================================================


# ------------------------------------------------------------
# FCM matching
# ------------------------------------------------------------

tab_fcm <- table(
  factor(
    labels_km,
    levels = seq_len(K)
  ),
  factor(
    labels_fcm,
    levels = seq_len(K)
  )
)

assignment_fcm <- clue::solve_LSAP(
  tab_fcm,
  maximum = TRUE
)

fcm_cluster_for_reference1 <- as.integer(
  assignment_fcm[1]
)


# ------------------------------------------------------------
# im.fuzzy matching
# ------------------------------------------------------------

tab_imf <- table(
  factor(
    labels_km,
    levels = seq_len(K)
  ),
  factor(
    labels_imf,
    levels = seq_len(K)
  )
)

assignment_imf <- clue::solve_LSAP(
  tab_imf,
  maximum = TRUE
)

imf_cluster_for_reference1 <- as.integer(
  assignment_imf[1]
)


cat(
  "\nFCM cluster corresponding to k-means Cluster 1:",
  fcm_cluster_for_reference1,
  "\n"
)

cat(
  "im.fuzzy cluster corresponding to k-means Cluster 1:",
  imf_cluster_for_reference1,
  "\n"
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
# 11. BUILD FCM MEMBERSHIP RASTER - CLUSTER 1
# ============================================================

r_fcm_membership1 <- template

vals_fcm <- rep(
  NA_real_,
  terra::ncell(template)
)

vals_fcm[valid_idx] <-
  U_fcm[
    ,
    fcm_cluster_for_reference1
  ]

terra::values(
  r_fcm_membership1
) <- vals_fcm


# ============================================================
# 12. BUILD im.fuzzy MEMBERSHIP RASTER - CLUSTER 1
# ============================================================

r_imf_membership1 <- template

vals_imf <- rep(
  NA_real_,
  terra::ncell(template)
)

vals_imf[valid_idx] <-
  U_imf[
    ,
    imf_cluster_for_reference1
  ]

terra::values(
  r_imf_membership1
) <- vals_imf


# ============================================================
# 13. ABSOLUTE MEMBERSHIP DIFFERENCE
# ============================================================

r_membership_difference <- abs(
  r_fcm_membership1 -
  r_imf_membership1
)


# ============================================================
# 14. NUMERICAL DIFFERENCE SUMMARY
# ============================================================

membership_difference <- abs(
  U_fcm[
    ,
    fcm_cluster_for_reference1
  ] -
  U_imf[
    ,
    imf_cluster_for_reference1
  ]
)


cat("\n")
cat("============================================\n")
cat("MEMBERSHIP DIFFERENCE SUMMARY\n")
cat("============================================\n")

cat(
  "Mean absolute difference:",
  mean(
    membership_difference,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median absolute difference:",
  median(
    membership_difference,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum absolute difference:",
  max(
    membership_difference,
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

difference_colors <- viridisLite::magma(
  100
)


# ============================================================
# 16. FINAL FIGURE
# ============================================================

png(
  filename = "monsummano_fuzzy_comparison_K2.png",
  width = 3000,
  height = 3600,
  res = 300
)

par(
  mfrow = c(3, 2),
  mar = c(1, 1, 3, 2)
)


# ------------------------------------------------------------
# (a) RGB image
# ------------------------------------------------------------

terra::plotRGB(
  img,
  r = 1,
  g = 2,
  b = 3,
  stretch = "lin",
  axes = FALSE,
  main = "(a) RGB image"
)


# ------------------------------------------------------------
# (b) k-means crisp partition
# ------------------------------------------------------------

terra::plot(
  r_km,
  col = cluster_colors,
  breaks = seq(
    0.5,
    K + 0.5,
    by = 1
  ),
  axes = FALSE,
  legend = FALSE,
  main = "(b) k-means partition"
)


# ------------------------------------------------------------
# (c) FCM membership
# ------------------------------------------------------------

terra::plot(
  r_fcm_membership1,
  col = membership_colors,
  range = c(0, 1),
  axes = FALSE,
  main = "(c) fuzzy C-means membership - Cluster 1"
)


# ------------------------------------------------------------
# (d) im.fuzzy membership
# ------------------------------------------------------------

terra::plot(
  r_imf_membership1,
  col = membership_colors,
  range = c(0, 1),
  axes = FALSE,
  main = "(d) im.fuzzy() membership - Cluster 1"
)


# ------------------------------------------------------------
# (e) Absolute membership difference
# ------------------------------------------------------------

terra::plot(
  r_membership_difference,
  col = difference_colors,
  axes = FALSE,
  main = "(e) absolute membership difference"
)


# ------------------------------------------------------------
# Empty sixth panel
# ------------------------------------------------------------

plot.new()


dev.off()


# ============================================================
# 17. SAVE RASTERS
# ============================================================

terra::writeRaster(
  r_km,
  "monsummano_kmeans_partition_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_fcm_membership1,
  "monsummano_fcm_membership_cluster1_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_imf_membership1,
  "monsummano_imfuzzy_membership_cluster1_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_membership_difference,
  "monsummano_membership_difference_K2.tif",
  overwrite = TRUE
)


# ============================================================
# END
# ============================================================
