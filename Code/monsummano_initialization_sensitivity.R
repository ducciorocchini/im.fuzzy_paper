# ============================================================
# MONSUMMANO - INITIALIZATION SENSITIVITY
#
# Purpose:
# Assess the sensitivity of im.fuzzy() to k-means initialization
# by repeating the analysis across multiple random seeds.
#
# Metrics:
#   - FPC: Fuzzy Partition Coefficient
#   - PE: Partition Entropy
#   - ARI: agreement with a reference partition
#   - NMI: agreement with a reference partition
#
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(mclust)
library(aricode)
library(ggplot2)


# ------------------------------------------------------------
# 2. Load Monsummano image
# ------------------------------------------------------------

img <- rast(
  "https://raw.githubusercontent.com/ducciorocchini/im.fuzzy_paper/main/Data/monsummano.png"
)

print(img)


# ------------------------------------------------------------
# 3. Analysis parameters
# ------------------------------------------------------------

K <- 4
m <- 2

# Number of random initializations
n_runs <- 30

# Seeds to test
seeds <- 1:n_runs


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
# 5. Helper functions
# ------------------------------------------------------------

# Fuzzy Partition Coefficient
FPC <- function(U) {

  mean(
    rowSums(U^2)
  )
}


# Partition Entropy
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
# 6. Reference partition
# ------------------------------------------------------------
#
# We use seed = 1 as a reference run.
#
# ARI and NMI measure agreement between each subsequent run
# and this reference partition.
#
# These are NOT accuracy metrics because no external ground
# truth exists for the Monsummano image.
# ------------------------------------------------------------

reference <- im.fuzzy(
  img,
  num_clusters = K,
  m = m,
  seed = 1,
  do_plot = FALSE
)

reference_membership <- terra::values(
  reference$memberships,
  mat = TRUE
)

reference_membership <- reference_membership[
  valid_idx,
  ,
  drop = FALSE
]

reference_labels <- apply(
  reference_membership,
  1,
  which.max
)


# ------------------------------------------------------------
# 7. Storage object
# ------------------------------------------------------------

results <- data.frame(
  Seed = integer(),
  FPC = numeric(),
  PE = numeric(),
  ARI = numeric(),
  NMI = numeric()
)


# ============================================================
# 8. Repeat im.fuzzy across seeds
# ============================================================

for (s in seeds) {

  cat(
    "Running seed:",
    s,
    "\n"
  )

  # ----------------------------------------------------------
  # Run im.fuzzy
  # ----------------------------------------------------------

  res <- im.fuzzy(
    img,
    num_clusters = K,
    m = m,
    seed = s,
    do_plot = FALSE
  )


  # ----------------------------------------------------------
  # Extract membership matrix
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
  # Crisp labels for agreement metrics
  # ----------------------------------------------------------

  labels <- apply(
    U,
    1,
    which.max
  )


  # ----------------------------------------------------------
  # Calculate metrics
  # ----------------------------------------------------------

  fpc_value <- FPC(U)

  pe_value <- partition_entropy(U)

  ari_value <- mclust::adjustedRandIndex(
    labels,
    reference_labels
  )

  nmi_value <- aricode::NMI(
    labels,
    reference_labels
  )


  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------

  results <- rbind(
    results,
    data.frame(
      Seed = s,
      FPC = fpc_value,
      PE = pe_value,
      ARI = ari_value,
      NMI = nmi_value
    )
  )
}


# ============================================================
# 9. Inspect results
# ============================================================

print(results)


# ------------------------------------------------------------
# 10. Summary statistics
# ------------------------------------------------------------

summary_results <- data.frame(

  Metric = c(
    "FPC",
    "Partition Entropy",
    "ARI",
    "NMI"
  ),

  Mean = c(
    mean(results$FPC),
    mean(results$PE),
    mean(results$ARI),
    mean(results$NMI)
  ),

  SD = c(
    sd(results$FPC),
    sd(results$PE),
    sd(results$ARI),
    sd(results$NMI)
  ),

  Min = c(
    min(results$FPC),
    min(results$PE),
    min(results$ARI),
    min(results$NMI)
  ),

  Max = c(
    max(results$FPC),
    max(results$PE),
    max(results$ARI),
    max(results$NMI)
  )
)

cat("\nSUMMARY OF INITIALIZATION SENSITIVITY\n")

print(summary_results)


# ============================================================
# 11. Prepare fuzzy metrics for plotting
# ============================================================

fuzzy_long <- rbind(

  data.frame(
    Seed = results$Seed,
    Metric = "FPC",
    Value = results$FPC
  ),

  data.frame(
    Seed = results$Seed,
    Metric = "Partition Entropy",
    Value = results$PE
  )
)


# ============================================================
# 12. GRAPH 1
# Fuzzy metrics across initializations
# ============================================================

p1 <- ggplot(
  fuzzy_long,
  aes(
    x = factor(Seed),
    y = Value,
    fill = Metric
  )
) +

  geom_col(
    position = "dodge"
  ) +

  scale_fill_manual(
    values = c(
      "#440154FF",
      "#FDE725FF"
    )
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Sensitivity of fuzzy partition characteristics",
    subtitle = "Monsummano image - 30 random initializations",
    x = "Random seed",
    y = "Index value"
  ) +

  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

print(p1)


# ============================================================
# 13. Prepare agreement metrics for plotting
# ============================================================

agreement_long <- rbind(

  data.frame(
    Seed = results$Seed,
    Metric = "ARI",
    Value = results$ARI
  ),

  data.frame(
    Seed = results$Seed,
    Metric = "NMI",
    Value = results$NMI
  )
)


# ============================================================
# 14. GRAPH 2
# Agreement with reference initialization
# ============================================================

p2 <- ggplot(
  agreement_long,
  aes(
    x = factor(Seed),
    y = Value,
    fill = Metric
  )
) +

  geom_col(
    position = "dodge"
  ) +

  scale_fill_manual(
    values = c(
      "#440154FF",
      "#FDE725FF"
    )
  ) +

  coord_cartesian(
    ylim = c(0, 1)
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Sensitivity to k-means initialization",
    subtitle = "Agreement with the reference partition (seed = 1)",
    x = "Random seed",
    y = "Agreement"
  ) +

  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

print(p2)


# ============================================================
# 15. GRAPH 3
# Boxplots summarizing stability
# ============================================================

stability_long <- rbind(

  data.frame(
    Metric = "ARI",
    Value = results$ARI
  ),

  data.frame(
    Metric = "NMI",
    Value = results$NMI
  ),

  data.frame(
    Metric = "FPC",
    Value = results$FPC
  ),

  data.frame(
    Metric = "PE",
    Value = results$PE
  )
)


p3 <- ggplot(
  stability_long,
  aes(
    x = Metric,
    y = Value,
    fill = Metric
  )
) +

  geom_boxplot(
    width = 0.6
  ) +

  scale_fill_viridis_d() +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Overall stability across random initializations",
    subtitle = "Monsummano image - 30 runs",
    x = NULL,
    y = "Metric value"
  ) +

  theme(
    legend.position = "none"
  )

print(p3)


# ============================================================
# 16. Save numerical results
# ============================================================

write.csv(
  results,
  "monsummano_initialization_results.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "monsummano_initialization_summary.csv",
  row.names = FALSE
)


# ============================================================
# 17. Save figures
# ============================================================

ggsave(
  "monsummano_initialization_fuzzy_metrics.png",
  p1,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_initialization_agreement.png",
  p2,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_initialization_stability.png",
  p3,
  width = 7,
  height = 5,
  dpi = 300
)

