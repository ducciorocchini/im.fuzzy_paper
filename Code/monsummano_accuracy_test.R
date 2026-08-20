# ============================================================
# MONSUMMANO VALIDATION EXPERIMENT
#
# Comparison among:
#   1. k-means
#   2. fuzzy c-means
#   3. im.fuzzy()
#
# Metrics:
#   - ARI and NMI: agreement among clustering solutions
#   - Silhouette and Dunn: internal cluster validity
#   - FPC and Partition Entropy: fuzzy partition quality
#
# IMPORTANT:
# Monsummano has no reference thematic classification.
# Therefore ARI and NMI are used as measures of agreement
# between methods, NOT as accuracy against ground truth.
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(e1071)     # fuzzy c-means
library(mclust)    # Adjusted Rand Index
library(aricode)   # Normalized Mutual Information
library(cluster)   # Silhouette
library(ggplot2)
library(imageRy)


# If needed:
# install.packages(c("terra", "e1071", "mclust",
#                    "aricode", "cluster", "ggplot2"))


# ------------------------------------------------------------
# 2. Load the Monsummano image
# ------------------------------------------------------------

img <- rast("https://raw.githubusercontent.com/ducciorocchini/im.fuzzy_paper/main/Data/monsummano.png") 
img <- flip(img)
img <- c(img[[1]], img[[2]], img[[3]])
plot(img)

im.plotRGB(img, 1, 2, 3)

# ------------------------------------------------------------
# 3. Extract pixel values
# ------------------------------------------------------------

# Rows = pixels
# Columns = raster bands
X <- terra::as.matrix(img)

# Retain only pixels having valid values in ALL bands
valid_idx <- complete.cases(X)

X_valid <- X[valid_idx, , drop = FALSE]

cat("Number of valid pixels:", nrow(X_valid), "\n")
cat("Number of bands:", ncol(X_valid), "\n")


# ------------------------------------------------------------
# 4. Scale every band independently to 0-255
# ------------------------------------------------------------
#
# This reproduces the preprocessing performed internally by
# im.fuzzy(), ensuring that all methods operate in the same
# feature space.
#

minmax_255 <- function(x) {

  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)

  # Constant bands contain no discriminatory information
  if (xmax == xmin) {
    return(rep(0, length(x)))
  }

  255 * (x - xmin) / (xmax - xmin)
}


X_scaled <- apply(
  X_valid,
  2,
  minmax_255
)

# Ensure matrix structure is retained
X_scaled <- as.matrix(X_scaled)


# ------------------------------------------------------------
# 5. Analysis parameters
# ------------------------------------------------------------

K <- 4
seed <- 42
m <- 2

set.seed(seed)


# ============================================================
# 6. K-MEANS
# ============================================================

km <- kmeans(
  X_scaled,
  centers = K,
  nstart = 25
)

km_labels <- km$cluster


# ============================================================
# 7. FUZZY C-MEANS
# ============================================================

set.seed(seed)

fcm <- e1071::cmeans(
  X_scaled,
  centers = K,
  m = m
)

# Fuzzy membership matrix
fcm_membership <- fcm$membership

# Crisp version obtained from maximum membership.
# This is needed only for ARI, NMI, Silhouette and Dunn.
fcm_labels <- apply(
  fcm_membership,
  1,
  which.max
)


# ============================================================
# 8. im.fuzzy()
# ============================================================

set.seed(seed)

res <- im.fuzzy(
  img,
  num_clusters = K,
  m = m,
  seed = seed,
  do_plot = FALSE
)

# Extract fuzzy membership values
imf_membership_all <- terra::values(
  res$memberships,
  mat = TRUE
)

# Retain the same valid pixels used by the other methods
imf_membership <- imf_membership_all[
  valid_idx,
  ,
  drop = FALSE
]

# Defuzzification by maximum membership.
# Again, this is ONLY for metrics requiring crisp partitions.
imf_labels <- apply(
  imf_membership,
  1,
  which.max
)


# ============================================================
# 9. PAIRWISE AGREEMENT
#
# ARI and NMI are label-invariant measures.
#
# They quantify similarity between cluster partitions.
#
# They are NOT interpreted here as classification accuracy
# because no independent ground truth exists for Monsummano.
# ============================================================

agreement <- data.frame(

  Comparison = c(
    "im.fuzzy vs k-means",
    "im.fuzzy vs fuzzy c-means",
    "fuzzy c-means vs k-means"
  ),

  ARI = c(

    mclust::adjustedRandIndex(
      imf_labels,
      km_labels
    ),

    mclust::adjustedRandIndex(
      imf_labels,
      fcm_labels
    ),

    mclust::adjustedRandIndex(
      fcm_labels,
      km_labels
    )
  ),

  NMI = c(

    aricode::NMI(
      imf_labels,
      km_labels
    ),

    aricode::NMI(
      imf_labels,
      fcm_labels
    ),

    aricode::NMI(
      fcm_labels,
      km_labels
    )
  )
)

cat("\nPAIRWISE AGREEMENT\n")
print(agreement)

#       Comparison       ARI       NMI
#1       im.fuzzy vs k-means 0.9930879 0.9868549
#2 im.fuzzy vs fuzzy c-means 0.9369748 0.9285782
#3  fuzzy c-means vs k-means 0.9343204 0.9264487

# ============================================================
# 10. FUZZY PARTITION COEFFICIENT
#
# FPC = (1/N) * sum_i sum_k u_ik^2
#
# Higher values indicate a sharper fuzzy partition.
#
# Minimum theoretical value = 1/K
# Maximum = 1
# ============================================================

FPC <- function(U) {

  mean(
    rowSums(U^2)
  )
}


# ============================================================
# 11. PARTITION ENTROPY
#
# PE = -(1/N) * sum_i sum_k u_ik log(u_ik)
#
# Lower values indicate a sharper partition.
# Higher values indicate greater ambiguity/fuzziness.
#
# Maximum theoretical value = log(K)
# ============================================================

partition_entropy <- function(U) {

  # Avoid log(0)
  U_safe <- pmax(U, .Machine$double.eps)

  -mean(
    rowSums(
      U_safe * log(U_safe)
    )
  )
}


# ------------------------------------------------------------
# 12. Calculate fuzzy quality metrics
# ------------------------------------------------------------

fuzzy_quality <- data.frame(

  Method = c(
    "fuzzy c-means",
    "im.fuzzy"
  ),

  FPC = c(
    FPC(fcm_membership),
    FPC(imf_membership)
  ),

  PE = c(
    partition_entropy(fcm_membership),
    partition_entropy(imf_membership)
  )
)

cat("\nFUZZY PARTITION QUALITY\n")
print(fuzzy_quality)

#         Method       FPC        PE
#1 fuzzy c-means 0.7323550 0.5044475
#2      im.fuzzy 0.7329669 0.5047415

# ============================================================
# 13. INTERNAL VALIDITY
#
# Silhouette and Dunn require pairwise distances.
#
# Computing a complete N x N distance matrix for a large image
# would be unnecessarily expensive.
#
# We therefore calculate these metrics on the SAME reproducible
# random sample of pixels for all methods.
# ============================================================

sample_size <- min(
  3000,
  nrow(X_scaled)
)

set.seed(seed)

sample_idx <- sample(
  seq_len(nrow(X_scaled)),
  size = sample_size,
  replace = FALSE
)

X_sample <- X_scaled[
  sample_idx,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 14. Distance matrix for sampled pixels
# ------------------------------------------------------------

D <- dist(X_sample)


# ------------------------------------------------------------
# 15. Mean Silhouette
#
# Higher values indicate better compactness/separation.
# ------------------------------------------------------------

mean_silhouette <- function(labels, D) {

  s <- cluster::silhouette(
    labels,
    D
  )

  mean(s[, "sil_width"])
}


# ------------------------------------------------------------
# 16. Dunn Index
#
# Dunn =
# minimum inter-cluster distance
# --------------------------------
# maximum intra-cluster diameter
#
# Higher values indicate better separated and more compact
# clusters.
# ------------------------------------------------------------

dunn_index <- function(X, labels) {

  D <- as.matrix(
    dist(X)
  )

  groups <- unique(labels)

  # --------------------------
  # Maximum within-cluster
  # diameter
  # --------------------------

  diameters <- sapply(
    groups,
    function(g) {

      idx <- which(labels == g)

      if (length(idx) <= 1) {
        return(0)
      }

      max(
        D[idx, idx, drop = FALSE]
      )
    }
  )

  max_diameter <- max(diameters)


  # --------------------------
  # Minimum distance between
  # two different clusters
  # --------------------------

  min_intercluster <- Inf

  for (i in seq_len(length(groups) - 1)) {

    for (j in (i + 1):length(groups)) {

      idx_i <- which(
        labels == groups[i]
      )

      idx_j <- which(
        labels == groups[j]
      )

      current_min <- min(
        D[idx_i, idx_j, drop = FALSE]
      )

      if (current_min < min_intercluster) {

        min_intercluster <- current_min

      }
    }
  }

  min_intercluster / max_diameter
}


# ------------------------------------------------------------
# 17. Calculate internal validity metrics
# ------------------------------------------------------------

internal_validity <- data.frame(

  Method = c(
    "k-means",
    "fuzzy c-means",
    "im.fuzzy"
  ),

  Silhouette = c(

    mean_silhouette(
      km_labels[sample_idx],
      D
    ),

    mean_silhouette(
      fcm_labels[sample_idx],
      D
    ),

    mean_silhouette(
      imf_labels[sample_idx],
      D
    )
  ),

  Dunn = c(

    dunn_index(
      X_sample,
      km_labels[sample_idx]
    ),

    dunn_index(
      X_sample,
      fcm_labels[sample_idx]
    ),

    dunn_index(
      X_sample,
      imf_labels[sample_idx]
    )
  )
)

cat("\nINTERNAL CLUSTER VALIDITY\n")
print(internal_validity)

#         Method Silhouette        Dunn
#1       k-means  0.4808493 0.006927926
#2 fuzzy c-means  0.4785145 0.006957711
#3      im.fuzzy  0.4806610 0.006858975

# ============================================================
# 18. GRAPH 1
# Pairwise ARI and NMI
# ============================================================

agreement_long <- rbind(

  data.frame(
    Comparison = agreement$Comparison,
    Metric = "ARI",
    Value = agreement$ARI
  ),

  data.frame(
    Comparison = agreement$Comparison,
    Metric = "NMI",
    Value = agreement$NMI
  )
)


p1 <- ggplot(
  agreement_long,
  aes(
    x = Comparison,
    y = Value,
    fill = Metric
  )
) +

  geom_col(
    position = "dodge"
  ) +

  coord_cartesian(
    ylim = c(0, 1)
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Agreement among clustering methods",
    subtitle = "Monsummano image",
    x = NULL,
    y = "Agreement"
  ) +

  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  ) +

  scale_fill_viridis_d(
    option = "viridis",
    end = 0.9
  )


print(p1)


# ============================================================
# 19. GRAPH 2
# Fuzzy partition quality
#
# FPC and PE have opposite interpretations:
#
# higher FPC = sharper partition
# lower PE   = sharper partition
#
# Therefore they should be interpreted separately.
# ============================================================

fuzzy_long <- rbind(

  data.frame(
    Method = fuzzy_quality$Method,
    Metric = "FPC",
    Value = fuzzy_quality$FPC
  ),

  data.frame(
    Method = fuzzy_quality$Method,
    Metric = "Partition Entropy",
    Value = fuzzy_quality$PE
  )
)


p2 <- ggplot(
  fuzzy_long,
  aes(
    x = Method,
    y = Value,
    fill = Metric
  )
) +

  geom_col(
    position = "dodge"
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Fuzzy partition characteristics",
    subtitle = "Monsummano image",
    x = NULL,
    y = "Index value"
  ) +

  scale_fill_viridis_d(
    option = "viridis",
    end = 0.9
  )

print(p2)


# ============================================================
# 20. GRAPH 3
# Internal cluster validity
# ============================================================

internal_long <- rbind(

  data.frame(
    Method = internal_validity$Method,
    Metric = "Silhouette",
    Value = internal_validity$Silhouette
  ),

  data.frame(
    Method = internal_validity$Method,
    Metric = "Dunn",
    Value = internal_validity$Dunn
  )
)


p3 <- ggplot(
  internal_long,
  aes(
    x = Method,
    y = Value,
    fill = Metric
  )
) +

  geom_col(
    position = "dodge"
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Internal cluster validity",
    subtitle = paste(
      "Monsummano image -",
      sample_size,
      "randomly sampled pixels"
    ),
    x = NULL,
    y = "Index value"
  ) +

  scale_fill_viridis_d(
    option = "viridis",
    end = 0.9
  )

print(p3)


# ============================================================
# 21. Save results
# ============================================================

write.csv(
  agreement,
  "monsummano_agreement.csv",
  row.names = FALSE
)

write.csv(
  fuzzy_quality,
  "monsummano_fuzzy_quality.csv",
  row.names = FALSE
)

write.csv(
  internal_validity,
  "monsummano_internal_validity.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Save figures
# ------------------------------------------------------------

ggsave(
  "monsummano_agreement.png",
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_fuzzy_quality.png",
  p2,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_internal_validity.png",
  p3,
  width = 7,
  height = 5,
  dpi = 300
)
