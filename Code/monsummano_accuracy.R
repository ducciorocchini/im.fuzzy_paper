# ============================================================
# Synthetic ground-truth validation using monsummano.png
# Comparison among k-means, fuzzy c-means, and im.fuzzy()
# ============================================================

# -----------------------------
# 1. Load required packages
# -----------------------------

library(terra)
library(ggplot2)
library(mclust)
library(aricode)
library(e1071)

# -----------------------------
# 2. Load the image
# -----------------------------

# Set your working directory if needed:
# setwd("~/path/to/your/folder")

img <- rast("monsummano.png")

# Optional: inspect image
print(img)
plotRGB(img, r = 1, g = 2, b = 3, stretch = "lin")

# -----------------------------
# 3. Extract pixel values
# -----------------------------

X <- as.matrix(img)

# Remove NA pixels
valid_idx <- complete.cases(X)
X_valid <- X[valid_idx, ]

# Scale each band to 0–255
X_scaled <- scale(X_valid)

# -----------------------------
# 4. Set analysis parameters
# -----------------------------

K <- 4
seed <- 42

set.seed(seed)

# -----------------------------
# 5. Create pseudo ground truth
# -----------------------------
# Since monsummano.png has no true labels, we create a reference
# classification using k-means. This is NOT real ground truth,
# but it allows method agreement to be assessed quantitatively.

ref <- kmeans(X_scaled, centers = K, nstart = 25)
true_labels <- ref$cluster

# -----------------------------
# 6. Method 1: k-means
# -----------------------------

km <- kmeans(X_scaled, centers = K, nstart = 25)
km_labels <- km$cluster

# -----------------------------
# 7. Method 2: fuzzy c-means
# -----------------------------

fcm <- cmeans(X_scaled, centers = K, m = 2)
fcm_labels <- apply(fcm$membership, 1, which.max)

# -----------------------------
# 8. Method 3: im.fuzzy()
# -----------------------------
# Requires imageRy development version:
# devtools::install_github("ducciorocchini/imageRy")

res <- im.fuzzy(img, num_clusters = K, seed = seed, do_plot = FALSE)

mem <- values(res$memberships)
mem_valid <- mem[valid_idx, ]

imf_labels <- apply(mem_valid, 1, which.max)

# -----------------------------
# 9. Compute validation metrics
# -----------------------------

results <- data.frame(
  Method = c("k-means", "fuzzy c-means", "im.fuzzy"),
  ARI = c(
    adjustedRandIndex(true_labels, km_labels),
    adjustedRandIndex(true_labels, fcm_labels),
    adjustedRandIndex(true_labels, imf_labels)
  ),
  NMI = c(
    NMI(true_labels, km_labels),
    NMI(true_labels, fcm_labels),
    NMI(true_labels, imf_labels)
  )
)

print(results)

# -----------------------------
# 10. Convert results to long format
# -----------------------------

results_long <- data.frame(
  Method = rep(results$Method, 2),
  Metric = rep(c("ARI", "NMI"), each = nrow(results)),
  Value = c(results$ARI, results$NMI)
)

# -----------------------------
# 11. Plot ARI and NMI
# -----------------------------

p <- ggplot(results_long, aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Clustering agreement on monsummano.png",
    subtitle = "Comparison using pseudo-reference labels",
    x = "Method",
    y = "Score"
  ) +
  ylim(0, 1)

print(p)

# -----------------------------
# 12. Save plot
# -----------------------------

# ggsave(
#  filename = "monsummano_validation_metrics.png",
#  plot = p,
#  width = 8,
#  height = 5,
#  dpi = 300
# )
