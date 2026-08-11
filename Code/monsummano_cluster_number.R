# ============================================================
# MONSUMMANO - SENSITIVITY TO NUMBER OF CLUSTERS
#
# Purpose:
# Evaluate how the choice of K affects:
#   - internal cluster validity
#   - fuzzy partition characteristics
#
# Metrics:
#   - Mean Silhouette
#   - Fuzzy Partition Coefficient (FPC)
#   - Partition Entropy (PE)
#
# Candidate number of clusters:
#   K = 2, ..., 10
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(cluster)
library(ggplot2)
library(imageRy)
library(dplyr)


# ------------------------------------------------------------
# 2. Load Monsummano image
# ------------------------------------------------------------

img <- rast(
  "https://raw.githubusercontent.com/ducciorocchini/im.fuzzy_paper/main/Data/monsummano.png"
)

print(img)


# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------

K_values <- 2:10

m <- 2
seed <- 42

# Silhouette requires a distance matrix.
# Use a fixed reproducible sample to avoid excessive memory use.
sample_size <- 3000


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
# 5. Scale bands to 0-255
#
# Same feature space used internally by im.fuzzy()
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


# ------------------------------------------------------------
# 6. Fixed random sample for Silhouette
# ------------------------------------------------------------

set.seed(seed)

sample_size <- min(
  sample_size,
  nrow(X_scaled)
)

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

D_sample <- dist(X_sample)


# ============================================================
# 7. Helper functions
# ============================================================

# ------------------------------------------------------------
# Fuzzy Partition Coefficient
#
# Higher values = sharper memberships
# ------------------------------------------------------------

FPC <- function(U) {

  mean(
    rowSums(U^2)
  )
}


# ------------------------------------------------------------
# Partition Entropy
#
# Lower values = sharper memberships
# ------------------------------------------------------------

partition_entropy <- function(U) {

  U_safe <- pmax(
    U,
    .Machine$double.eps
  )

  -mean(
    rowSums(
      U_safe * log(U_safe)
    )
  )
}


# ------------------------------------------------------------
# Mean Silhouette
# ------------------------------------------------------------

mean_silhouette <- function(labels, D) {

  sil <- cluster::silhouette(
    labels,
    D
  )

  mean(
    sil[, "sil_width"]
  )
}


# ============================================================
# 8. Storage object
# ============================================================

results <- data.frame(
  K = integer(),
  Silhouette = numeric(),
  FPC = numeric(),
  PE = numeric()
)


# ============================================================
# 9. Loop over candidate K values
# ============================================================

for (K in K_values) {

  cat(
    "Running K =",
    K,
    "\n"
  )


  # ----------------------------------------------------------
  # Run im.fuzzy
  # ----------------------------------------------------------

  res <- imageRy::im.fuzzy(
    img,
    num_clusters = K,
    m = m,
    seed = seed,
    do_plot = FALSE
  )


  # ----------------------------------------------------------
  # Extract memberships
  # ----------------------------------------------------------

  U_all <- terra::values(
    res$memberships,
    mat = TRUE
  )

  U <- U_all[
    valid_idx,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Defuzzify using maximum membership
  # ----------------------------------------------------------

  labels <- apply(
    U,
    1,
    which.max
  )


  # ----------------------------------------------------------
  # Calculate fuzzy metrics
  # ----------------------------------------------------------

  fpc_value <- FPC(U)

  pe_value <- partition_entropy(U)


  # ----------------------------------------------------------
  # Calculate Silhouette on fixed sample
  # ----------------------------------------------------------

  silhouette_value <- mean_silhouette(
    labels[sample_idx],
    D_sample
  )


  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------

  results <- rbind(
    results,
    data.frame(
      K = K,
      Silhouette = silhouette_value,
      FPC = fpc_value,
      PE = pe_value
    )
  )
}


# ============================================================
# 10. Inspect results
# ============================================================

cat("\n")
cat("============================================\n")
cat("RESULTS\n")
cat("============================================\n")

print(results)


# ============================================================
# 11. Identify best K according to each metric
# ============================================================

best_silhouette <- results$K[
  which.max(results$Silhouette)
]

best_fpc <- results$K[
  which.max(results$FPC)
]

best_pe <- results$K[
  which.min(results$PE)
]


cat("\n")
cat("============================================\n")
cat("BEST K BY METRIC\n")
cat("============================================\n")

cat(
  "Maximum Silhouette: K =",
  best_silhouette,
  "\n"
)

cat(
  "Maximum FPC: K =",
  best_fpc,
  "\n"
)

cat(
  "Minimum PE: K =",
  best_pe,
  "\n"
)


# ============================================================
# 12. Normalize FPC for comparison across different K
#
# Raw FPC is influenced by K because its theoretical minimum
# is 1/K.
#
# Normalized FPC:
#
#   FPC_norm = (FPC - 1/K) / (1 - 1/K)
#
# ranges approximately from 0 to 1.
# ============================================================

results <- results %>%
  mutate(
    FPC_normalized =
      (FPC - 1 / K) /
      (1 - 1 / K),

    PE_normalized =
      PE / log(K)
  )


print(results)


# ============================================================
# 13. GRAPH 1
# Mean Silhouette vs K
# ============================================================

p_silhouette <- ggplot(
  results,
  aes(
    x = K,
    y = Silhouette
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks = K_values
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Internal cluster validity across candidate K",
    subtitle = "Monsummano image",
    x = "Number of clusters (K)",
    y = "Mean Silhouette"
  )

print(p_silhouette)


# ============================================================
# 14. GRAPH 2
# Normalized fuzzy partition metrics vs K
# ============================================================

fuzzy_long <- rbind(

  data.frame(
    K = results$K,
    Metric = "Normalized FPC",
    Value = results$FPC_normalized
  ),

  data.frame(
    K = results$K,
    Metric = "Normalized PE",
    Value = results$PE_normalized
  )
)


p_fuzzy <- ggplot(
  fuzzy_long,
  aes(
    x = K,
    y = Value,
    color = Metric,
    group = Metric
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks = K_values
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Fuzzy partition characteristics across candidate K",
    subtitle = "Monsummano image",
    x = "Number of clusters (K)",
    y = "Normalized index value",
    color = "Metric"
  )

print(p_fuzzy)


# ============================================================
# 15. GRAPH 3
# Raw FPC and PE vs K
# ============================================================

raw_long <- rbind(

  data.frame(
    K = results$K,
    Metric = "FPC",
    Value = results$FPC
  ),

  data.frame(
    K = results$K,
    Metric = "Partition Entropy",
    Value = results$PE
  )
)


p_raw <- ggplot(
  raw_long,
  aes(
    x = K,
    y = Value,
    color = Metric,
    group = Metric
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks = K_values
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Raw fuzzy validity indices across candidate K",
    subtitle = "Monsummano image",
    x = "Number of clusters (K)",
    y = "Index value",
    color = "Metric"
  )

print(p_raw)


# ============================================================
# 16. Save numerical results
# ============================================================

write.csv(
  results,
  "monsummano_cluster_number_results.csv",
  row.names = FALSE
)


# ============================================================
# 17. Save figures
# ============================================================

ggsave(
  filename = "monsummano_K_silhouette.png",
  plot = p_silhouette,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "monsummano_K_fuzzy_normalized.png",
  plot = p_fuzzy,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "monsummano_K_fuzzy_raw.png",
  plot = p_raw,
  width = 8,
  height = 5,
  dpi = 300
)

