# ============================================================
# MONSUMMANO - VISUAL COMPARISON OF CLUSTERING PARTITIONS
#
# Panels:
# (a) Original RGB image
# (b) k-means partition
# (c) fuzzy C-means maximum-membership partition
# (d) im.fuzzy() maximum-membership partition
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

img <- terra::rast(monsummano_url)
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

X <- terra::as.matrix(img)

valid_idx <- complete.cases(X)

X_valid <- X[
  valid_idx,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Scale bands independently to 0-255
#
# Same preprocessing used internally by im.fuzzy()
# ------------------------------------------------------------

scale_255 <- function(x) {

  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)

  if (xmax == xmin) {
    return(rep(0, length(x)))
  }

  255 * (x - xmin) / (xmax - xmin)
}


X_scaled <- apply(
  X_valid,
  2,
  scale_255
)

X_scaled <- as.matrix(X_scaled)


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

labels_fcm <- apply(
  fcm$membership,
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
# 9. MATCH CLUSTER LABELS TO K-MEANS
#
# Cluster IDs are arbitrary. This step relabels fuzzy C-means
# and im.fuzzy() so that corresponding groups use the same
# colors in all panels.
# ============================================================

match_labels <- function(reference, target, K) {

  tab <- table(
    factor(reference, levels = seq_len(K)),
    factor(target, levels = seq_len(K))
  )

  assignment <- clue::solve_LSAP(
    tab,
    maximum = TRUE
  )

  mapping <- integer(K)

  for (ref_cluster in seq_len(K)) {

    target_cluster <- as.integer(
      assignment[ref_cluster]
    )

    mapping[target_cluster] <- ref_cluster
  }

  mapping[target]
}


labels_fcm_matched <- match_labels(
  reference = labels_km,
  target = labels_fcm,
  K = K
)

labels_imf_matched <- match_labels(
  reference = labels_km,
  target = labels_imf,
  K = K
)


# ============================================================
# 10. CHECK AGREEMENT AFTER LABEL MATCHING
# ============================================================

cat(
  "\nProportion identical: k-means vs fuzzy C-means =",
  mean(labels_km == labels_fcm_matched),
  "\n"
)

cat(
  "Proportion identical: k-means vs im.fuzzy =",
  mean(labels_km == labels_imf_matched),
  "\n"
)


# ============================================================
# 11. BUILD PARTITION RASTERS
# ============================================================

template <- img[[1]]


make_partition <- function(labels) {

  r <- template

  vals <- rep(
    NA_integer_,
    terra::ncell(template)
  )

  vals[valid_idx] <- labels

  terra::values(r) <- vals

  r
}


r_km <- make_partition(
  labels_km
)

r_fcm <- make_partition(
  labels_fcm_matched
)

r_imf <- make_partition(
  labels_imf_matched
)


# ============================================================
# 12. COLORBLIND-FRIENDLY CLUSTER PALETTE
#
# Same visual language used elsewhere in the paper
# ============================================================

cluster_colors <- viridisLite::viridis(
  K,
  option = "D",
  end = 0.9
)


# ============================================================
# 13. PLOT FIGURE
# ============================================================

png(
  filename = "monsummano_partitions_K2.png",
  width = 3000,
  height = 2400,
  res = 300
)

par(
  mfrow = c(2, 2),
  mar = c(1, 1, 3, 1)
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
# (b) k-means
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
  main = "(b) k-means"
)


# ------------------------------------------------------------
# (c) fuzzy C-means
# ------------------------------------------------------------

terra::plot(
  r_fcm,
  col = cluster_colors,
  breaks = seq(
    0.5,
    K + 0.5,
    by = 1
  ),
  axes = FALSE,
  legend = FALSE,
  main = "(c) fuzzy C-means"
)


# ------------------------------------------------------------
# (d) im.fuzzy()
# ------------------------------------------------------------

terra::plot(
  r_imf,
  col = cluster_colors,
  breaks = seq(
    0.5,
    K + 0.5,
    by = 1
  ),
  axes = FALSE,
  legend = FALSE,
  main = "(d) im.fuzzy()"
)


dev.off()


# ============================================================
# 14. OPTIONAL: SAVE PARTITION RASTERS
# ============================================================

terra::writeRaster(
  r_km,
  "monsummano_kmeans_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_fcm,
  "monsummano_fuzzy_cmeans_K2.tif",
  overwrite = TRUE
)

terra::writeRaster(
  r_imf,
  "monsummano_imfuzzy_K2.tif",
  overwrite = TRUE
)

