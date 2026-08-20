# ============================================================
# MONSUMMANO - VISUAL COMPARISON OF CLUSTERING PARTITIONS
#
# Panels:
# (a) Original RGB image
# (b) k-means
# (c) fuzzy c-means
# (d) im.fuzzy()
#
# K = 4
# ============================================================

library(terra)
library(e1071)
library(imageRy)
library(clue)

# ------------------------------------------------------------
# 1. Load Monsummano image
# ------------------------------------------------------------

# CHANGE THIS to the same import command/path already used
# in your Monsummano analyses

img <- rast("monsummano.tif")

print(img)

# ------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------

K <- 4
m <- 2
seed <- 42

# ------------------------------------------------------------
# 3. Extract valid pixels
# ------------------------------------------------------------

X <- terra::as.matrix(img)

valid_idx <- complete.cases(X)

X_valid <- X[
  valid_idx,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 4. Same 0-255 scaling used by im.fuzzy()
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
# 5. K-MEANS
# ============================================================

set.seed(seed)

km <- kmeans(
  X_scaled,
  centers = K
)

labels_km <- km$cluster

# ============================================================
# 6. FUZZY C-MEANS
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
# 7. im.fuzzy()
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
# 8. MATCH CLUSTER LABELS TO K-MEANS
#
# Cluster numbers are arbitrary.
# We relabel FCM and im.fuzzy so that corresponding spatial
# partitions use the same colors.
# ============================================================

match_labels <- function(reference, target, K) {

  tab <- table(
    factor(reference, levels = seq_len(K)),
    factor(target, levels = seq_len(K))
  )

  # Maximize overlap between target and reference clusters
  assignment <- clue::solve_LSAP(
    tab,
    maximum = TRUE
  )

  # assignment[i] = target cluster corresponding to
  # reference cluster i
  mapping <- integer(K)

  for (ref_cluster in seq_len(K)) {

    target_cluster <- assignment[ref_cluster]

    mapping[target_cluster] <- ref_cluster
  }

  mapping[target]
}

labels_fcm_matched <- match_labels(
  labels_km,
  labels_fcm,
  K
)

labels_imf_matched <- match_labels(
  labels_km,
  labels_imf,
  K
)

# ============================================================
# 9. BUILD PARTITION RASTERS
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
# 10. COLORBLIND-FRIENDLY CLUSTER PALETTE
# ============================================================

cluster_colors <- viridisLite::viridis(
  K,
  option = "D",
  end = 0.9
)

# ============================================================
# 11. PLOT
# ============================================================

png(
  "monsummano_partitions.png",
  width = 3000,
  height = 2600,
  res = 300
)

par(
  mfrow = c(2, 2),
  mar = c(1, 1, 3, 1)
)

# ------------------------------------------------------------
# (a) Original image
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
# (c) fuzzy c-means
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
  main = "(c) fuzzy c-means"
)

# ------------------------------------------------------------
# (d) im.fuzzy
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
# 12. OPTIONAL CHECK OF AGREEMENT AFTER LABEL MATCHING
# ============================================================

cat(
  "\nProportion identical: k-means vs fuzzy c-means =",
  mean(labels_km == labels_fcm_matched),
  "\n"
)

cat(
  "Proportion identical: k-means vs im.fuzzy =",
  mean(labels_km == labels_imf_matched),
  "\n"
)

# ============================================================
# END
# ============================================================
