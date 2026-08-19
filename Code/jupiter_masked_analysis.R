# ============================================================
# JUPITER ANALYSIS WITH BACKGROUND MASK
#
# Workflow:
# 1. Import Jupiter image
# 2. Mask dark background outside planetary disk
# 3. Evaluate K = 2,...,10
# 4. Select K using Silhouette + normalized FPC + normalized PE
# 5. Run im.fuzzy() on masked Jupiter image
# ============================================================

library(terra)
library(cluster)
library(ggplot2)
library(imageRy)
library(dplyr)


# ------------------------------------------------------------
# 1. Load Jupiter image
# ------------------------------------------------------------

jupiter_url <- paste0(
  "https://raw.githubusercontent.com/",
  "ducciorocchini/im.fuzzy_paper/main/Data/",
  "STScI-01EVSV9JYWYARC4A8336KMKNTQ.jpeg"
)

img <- terra::rast(jupiter_url)

# Same orientation used in the paper
img <- terra::flip(img)

print(img)


# ------------------------------------------------------------
# 2. Compute RGB brightness
# ------------------------------------------------------------

brightness <- (
  img[[1]] +
  img[[2]] +
  img[[3]]
) / 3

names(brightness) <- "brightness"


# ------------------------------------------------------------
# 3. Create planetary-disk mask
#
# Pixels darker than threshold are treated as background.
# ------------------------------------------------------------

background_threshold <- 10

planet_mask <- brightness > background_threshold


# ------------------------------------------------------------
# 4. Apply mask to all RGB bands
# ------------------------------------------------------------

img_masked <- terra::mask(
  img,
  planet_mask,
  maskvalues = 0
)

print(img_masked)


# ------------------------------------------------------------
# 5. Inspect mask
# ------------------------------------------------------------

terra::plotRGB(
  img_masked,
  r = 1,
  g = 2,
  b = 3,
  stretch = "lin"
)

terra::plot(
  planet_mask,
  main = "Jupiter planetary-disk mask"
)


# ============================================================
# 6. Candidate K analysis
# ============================================================

K_values <- 2:10
m <- 2
seed <- 42
sample_size <- 3000


# ------------------------------------------------------------
# Extract valid pixels
# ------------------------------------------------------------

X <- terra::as.matrix(img_masked)

valid_idx <- complete.cases(X)

X_valid <- X[
  valid_idx,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Scale to 0-255, as in im.fuzzy()
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
# Fixed sample for Silhouette
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


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

FPC <- function(U) {
  mean(rowSums(U^2))
}


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
# 7. Loop over K
# ============================================================

results <- data.frame(
  K = integer(),
  Silhouette = numeric(),
  FPC = numeric(),
  PE = numeric()
)


for (K in K_values) {

  cat(
    "\nRunning masked Jupiter: K =",
    K,
    "\n"
  )

  res <- imageRy::im.fuzzy(
    img_masked,
    num_clusters = K,
    m = m,
    seed = seed,
    do_plot = FALSE
  )


  U_all <- terra::values(
    res$memberships,
    mat = TRUE
  )

  U <- U_all[
    valid_idx,
    ,
    drop = FALSE
  ]


  labels <- apply(
    U,
    1,
    which.max
  )


  silhouette_value <- mean_silhouette(
    labels[sample_idx],
    D_sample
  )

  fpc_value <- FPC(U)

  pe_value <- partition_entropy(U)


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
# 8. Normalize fuzzy indices
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
# 9. Best K according to each metric
# ============================================================

cat(
  "\nBest Silhouette K =",
  results$K[which.max(results$Silhouette)],
  "\n"
)

cat(
  "Best normalized FPC K =",
  results$K[which.max(results$FPC_normalized)],
  "\n"
)

cat(
  "Best normalized PE K =",
  results$K[which.min(results$PE_normalized)],
  "\n"
)


# ============================================================
# 10. Plot Silhouette
# ============================================================

p_silhouette <- ggplot(
  results,
  aes(
    x = K,
    y = Silhouette
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = K_values
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Internal cluster validity across candidate K",
    subtitle = "Masked Jupiter image",
    x = "Number of clusters (K)",
    y = "Mean Silhouette"
  )

print(p_silhouette)


# ============================================================
# 11. Plot normalized fuzzy metrics
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
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = K_values
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Fuzzy partition characteristics across candidate K",
    subtitle = "Masked Jupiter image",
    x = "Number of clusters (K)",
    y = "Normalized index value",
    color = "Metric"
  )

print(p_fuzzy)


# ============================================================
# 12. Save results
# ============================================================

write.csv(
  results,
  "jupiter_masked_cluster_number_results.csv",
  row.names = FALSE
)

ggsave(
  "jupiter_masked_K_silhouette.png",
  p_silhouette,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "jupiter_masked_K_fuzzy_normalized.png",
  p_fuzzy,
  width = 8,
  height = 5,
  dpi = 300
)
