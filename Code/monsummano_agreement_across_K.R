# ============================================================
# MONSUMMANO - AGREEMENT AMONG METHODS ACROSS K
#
# Purpose:
# Test whether the agreement among:
#
#   1. k-means
#   2. fuzzy C-means
#   3. im.fuzzy()
#
# remains high as the number of clusters increases.
#
# Candidate K:
#   2,...,10
#
# Metrics:
#   - ARI
#   - NMI
#   - FPC
#   - Partition Entropy
#
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(e1071)
library(imageRy)
library(mclust)
library(aricode)
library(ggplot2)
library(dplyr)
library(viridisLite)


# ------------------------------------------------------------
# 2. Load Monsummano image
#
# IMPORTANT:
# use exactly the same image used in the other Monsummano tests
# ------------------------------------------------------------

monsummano_url <- paste0(
  "https://raw.githubusercontent.com/",
  "ducciorocchini/im.fuzzy_paper/main/Data/",
  "monsummano.png"
)

img <- terra::rast(monsummano_url)

print(img)


# ------------------------------------------------------------
# 3. Parameters
# ------------------------------------------------------------

K_values <- 2:10

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
# 5. Scale each band to 0-255
#
# Same preprocessing used by im.fuzzy()
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
      rep(0, length(x))
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
# 6. Helper functions
# ============================================================


# ------------------------------------------------------------
# Fuzzy Partition Coefficient
# ------------------------------------------------------------

FPC <- function(U) {

  mean(
    rowSums(U^2)
  )
}


# ------------------------------------------------------------
# Partition Entropy
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


# ============================================================
# 7. Storage objects
# ============================================================

agreement_results <- data.frame(

  K = integer(),

  ARI_imf_km = numeric(),
  NMI_imf_km = numeric(),

  ARI_imf_fcm = numeric(),
  NMI_imf_fcm = numeric(),

  ARI_fcm_km = numeric(),
  NMI_fcm_km = numeric()
)


fuzzy_results <- data.frame(

  K = integer(),

  FPC_imf = numeric(),
  FPC_fcm = numeric(),

  PE_imf = numeric(),
  PE_fcm = numeric()
)


# ============================================================
# 8. LOOP ACROSS K
# ============================================================

for (K in K_values) {

  cat(
    "\n============================================\n"
  )

  cat(
    "Running K =",
    K,
    "\n"
  )

  cat(
    "============================================\n"
  )


  # ==========================================================
  # 8.1 K-MEANS
  # ==========================================================

  set.seed(seed)

  km <- stats::kmeans(
    X_scaled,
    centers = K
  )

  labels_km <- km$cluster


  # ==========================================================
  # 8.2 FUZZY C-MEANS
  # ==========================================================

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


  # ==========================================================
  # 8.3 im.fuzzy()
  # ==========================================================

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


  # ==========================================================
  # 8.4 ARI
  # ==========================================================

  ARI_imf_km <- mclust::adjustedRandIndex(
    labels_imf,
    labels_km
  )

  ARI_imf_fcm <- mclust::adjustedRandIndex(
    labels_imf,
    labels_fcm
  )

  ARI_fcm_km <- mclust::adjustedRandIndex(
    labels_fcm,
    labels_km
  )


  # ==========================================================
  # 8.5 NMI
  # ==========================================================

  NMI_imf_km <- aricode::NMI(
    labels_imf,
    labels_km
  )

  NMI_imf_fcm <- aricode::NMI(
    labels_imf,
    labels_fcm
  )

  NMI_fcm_km <- aricode::NMI(
    labels_fcm,
    labels_km
  )


  # ==========================================================
  # 8.6 FUZZY METRICS
  # ==========================================================

  FPC_imf <- FPC(
    U_imf
  )

  FPC_fcm <- FPC(
    U_fcm
  )


  PE_imf <- partition_entropy(
    U_imf
  )

  PE_fcm <- partition_entropy(
    U_fcm
  )


  # ==========================================================
  # 8.7 STORE AGREEMENT
  # ==========================================================

  agreement_results <- rbind(

    agreement_results,

    data.frame(

      K = K,

      ARI_imf_km = ARI_imf_km,
      NMI_imf_km = NMI_imf_km,

      ARI_imf_fcm = ARI_imf_fcm,
      NMI_imf_fcm = NMI_imf_fcm,

      ARI_fcm_km = ARI_fcm_km,
      NMI_fcm_km = NMI_fcm_km
    )
  )


  # ==========================================================
  # 8.8 STORE FUZZY METRICS
  # ==========================================================

  fuzzy_results <- rbind(

    fuzzy_results,

    data.frame(

      K = K,

      FPC_imf = FPC_imf,
      FPC_fcm = FPC_fcm,

      PE_imf = PE_imf,
      PE_fcm = PE_fcm
    )
  )


  # ----------------------------------------------------------
  # Console summary
  # ----------------------------------------------------------

  cat(
    "ARI im.fuzzy vs k-means =",
    round(ARI_imf_km, 4),
    "\n"
  )

  cat(
    "ARI im.fuzzy vs FCM =",
    round(ARI_imf_fcm, 4),
    "\n"
  )

  cat(
    "NMI im.fuzzy vs k-means =",
    round(NMI_imf_km, 4),
    "\n"
  )

  cat(
    "NMI im.fuzzy vs FCM =",
    round(NMI_imf_fcm, 4),
    "\n"
  )


  # ----------------------------------------------------------
  # Clean objects before next K
  # ----------------------------------------------------------

  rm(
    km,
    fcm,
    imf,
    U_fcm,
    U_imf,
    U_imf_all
  )

  gc()
}


# ============================================================
# 9. PRINT RESULTS
# ============================================================

cat("\n")
cat("============================================\n")
cat("AGREEMENT ACROSS K\n")
cat("============================================\n")

print(
  agreement_results
)


cat("\n")
cat("============================================\n")
cat("FUZZY METRICS ACROSS K\n")
cat("============================================\n")

print(
  fuzzy_results
)


# ============================================================
# 10. LONG FORMAT - ARI
# ============================================================

ARI_long <- rbind(

  data.frame(
    K = agreement_results$K,
    Comparison = "im.fuzzy vs k-means",
    Value = agreement_results$ARI_imf_km
  ),

  data.frame(
    K = agreement_results$K,
    Comparison = "im.fuzzy vs fuzzy C-means",
    Value = agreement_results$ARI_imf_fcm
  ),

  data.frame(
    K = agreement_results$K,
    Comparison = "fuzzy C-means vs k-means",
    Value = agreement_results$ARI_fcm_km
  )
)


# ============================================================
# 11. LONG FORMAT - NMI
# ============================================================

NMI_long <- rbind(

  data.frame(
    K = agreement_results$K,
    Comparison = "im.fuzzy vs k-means",
    Value = agreement_results$NMI_imf_km
  ),

  data.frame(
    K = agreement_results$K,
    Comparison = "im.fuzzy vs fuzzy C-means",
    Value = agreement_results$NMI_imf_fcm
  ),

  data.frame(
    K = agreement_results$K,
    Comparison = "fuzzy C-means vs k-means",
    Value = agreement_results$NMI_fcm_km
  )
)


# ============================================================
# 12. PAPER COLORS
# ============================================================

comparison_colors <- viridisLite::viridis(
  3,
  option = "D",
  end = 0.9
)


# ============================================================
# 13. ARI ACROSS K
# ============================================================

p_ari <- ggplot(
  ARI_long,
  aes(
    x = K,
    y = Value,
    color = Comparison,
    group = Comparison
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  scale_color_manual(
    values = comparison_colors
  ) +

  scale_x_continuous(
    breaks = K_values
  ) +

  coord_cartesian(
    ylim = c(0, 1)
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Agreement among clustering methods across K",
    subtitle = "Adjusted Rand Index - Monsummano image",
    x = "Number of clusters (K)",
    y = "ARI",
    color = "Comparison"
  )


print(
  p_ari
)


# ============================================================
# 14. NMI ACROSS K
# ============================================================

p_nmi <- ggplot(
  NMI_long,
  aes(
    x = K,
    y = Value,
    color = Comparison,
    group = Comparison
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  scale_color_manual(
    values = comparison_colors
  ) +

  scale_x_continuous(
    breaks = K_values
  ) +

  coord_cartesian(
    ylim = c(0, 1)
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Agreement among clustering methods across K",
    subtitle = "Normalized Mutual Information - Monsummano image",
    x = "Number of clusters (K)",
    y = "NMI",
    color = "Comparison"
  )


print(
  p_nmi
)


# ============================================================
# 15. FPC COMPARISON ACROSS K
# ============================================================

FPC_long <- rbind(

  data.frame(
    K = fuzzy_results$K,
    Method = "im.fuzzy",
    Value = fuzzy_results$FPC_imf
  ),

  data.frame(
    K = fuzzy_results$K,
    Method = "fuzzy C-means",
    Value = fuzzy_results$FPC_fcm
  )
)


paper_colors <- viridisLite::viridis(
  2,
  option = "D",
  end = 0.9
)


p_fpc <- ggplot(
  FPC_long,
  aes(
    x = K,
    y = Value,
    color = Method,
    group = Method
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  scale_color_manual(
    values = paper_colors
  ) +

  scale_x_continuous(
    breaks = K_values
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Fuzzy Partition Coefficient across K",
    subtitle = "Monsummano image",
    x = "Number of clusters (K)",
    y = "FPC",
    color = "Method"
  )


print(
  p_fpc
)


# ============================================================
# 16. PARTITION ENTROPY ACROSS K
# ============================================================

PE_long <- rbind(

  data.frame(
    K = fuzzy_results$K,
    Method = "im.fuzzy",
    Value = fuzzy_results$PE_imf
  ),

  data.frame(
    K = fuzzy_results$K,
    Method = "fuzzy C-means",
    Value = fuzzy_results$PE_fcm
  )
)


p_pe <- ggplot(
  PE_long,
  aes(
    x = K,
    y = Value,
    color = Method,
    group = Method
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  scale_color_manual(
    values = paper_colors
  ) +

  scale_x_continuous(
    breaks = K_values
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Partition Entropy across K",
    subtitle = "Monsummano image",
    x = "Number of clusters (K)",
    y = "Partition Entropy",
    color = "Method"
  )


print(
  p_pe
)


# ============================================================
# 17. SAVE NUMERICAL RESULTS
# ============================================================

write.csv(
  agreement_results,
  "monsummano_agreement_across_K.csv",
  row.names = FALSE
)

write.csv(
  fuzzy_results,
  "monsummano_fuzzy_metrics_across_K.csv",
  row.names = FALSE
)


# ============================================================
# 18. SAVE FIGURES
# ============================================================

ggsave(
  "monsummano_ARI_across_K.png",
  p_ari,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_NMI_across_K.png",
  p_nmi,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_FPC_across_K.png",
  p_fpc,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "monsummano_PE_across_K.png",
  p_pe,
  width = 8,
  height = 5,
  dpi = 300
)

###########

# ============================================================
# PUBLICATION FIGURE
# AGREEMENT ACROSS K: ARI + NMI
#
# Single ggplot with facets
# ============================================================

library(ggplot2)
library(viridisLite)


# ------------------------------------------------------------
# 1. Build single long data frame
# ------------------------------------------------------------

agreement_plot <- rbind(

  data.frame(
    K = agreement_results$K,
    Metric = "Adjusted Rand Index",
    Comparison = "fuzzy C-means vs k-means",
    Value = agreement_results$ARI_fcm_km
  ),

  data.frame(
    K = agreement_results$K,
    Metric = "Adjusted Rand Index",
    Comparison = "im.fuzzy vs fuzzy C-means",
    Value = agreement_results$ARI_imf_fcm
  ),

  data.frame(
    K = agreement_results$K,
    Metric = "Adjusted Rand Index",
    Comparison = "im.fuzzy vs k-means",
    Value = agreement_results$ARI_imf_km
  ),

  data.frame(
    K = agreement_results$K,
    Metric = "Normalized Mutual Information",
    Comparison = "fuzzy C-means vs k-means",
    Value = agreement_results$NMI_fcm_km
  ),

  data.frame(
    K = agreement_results$K,
    Metric = "Normalized Mutual Information",
    Comparison = "im.fuzzy vs fuzzy C-means",
    Value = agreement_results$NMI_imf_fcm
  ),

  data.frame(
    K = agreement_results$K,
    Metric = "Normalized Mutual Information",
    Comparison = "im.fuzzy vs k-means",
    Value = agreement_results$NMI_imf_km
  )
)


# ------------------------------------------------------------
# 2. Factor order
# ------------------------------------------------------------

agreement_plot$Comparison <- factor(
  agreement_plot$Comparison,
  levels = c(
    "fuzzy C-means vs k-means",
    "im.fuzzy vs fuzzy C-means",
    "im.fuzzy vs k-means"
  )
)


agreement_plot$Metric <- factor(
  agreement_plot$Metric,
  levels = c(
    "Adjusted Rand Index",
    "Normalized Mutual Information"
  ),
  labels = c(
    "(a) Adjusted Rand Index",
    "(b) Normalized Mutual Information"
  )
)


# ------------------------------------------------------------
# 3. Colorblind-friendly colors
# ------------------------------------------------------------

comparison_colors <- viridisLite::viridis(
  3,
  option = "D",
  end = 0.9
)


# ============================================================
# 4. FINAL FIGURE
# ============================================================

p_agreement_K <- ggplot(
  agreement_plot,
  aes(
    x = K,
    y = Value,
    color = Comparison,
    group = Comparison
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  facet_wrap(
    ~ Metric,
    ncol = 1
  ) +

  scale_color_manual(
    values = comparison_colors
  ) +

  scale_x_continuous(
    breaks = 2:10
  ) +

  scale_y_continuous(
    breaks = seq(
      0,
      1,
      0.25
    )
  ) +

  coord_cartesian(
    ylim = c(0, 1)
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Agreement among clustering methods across K",
    subtitle = "Monsummano image",
    x = "Number of clusters (K)",
    y = "Agreement",
    color = "Comparison"
  ) +

  theme(
    legend.position = "bottom",

    strip.text = element_text(
      size = 14,
      face = "bold"
    ),

    panel.spacing = grid::unit(
      1.2,
      "lines"
    )
  )


# ============================================================
# 5. DISPLAY
# ============================================================

print(
  p_agreement_K
)


# ============================================================
# 6. SAVE
# ============================================================

ggsave(
  filename = "monsummano_agreement_across_K_combined.png",
  plot = p_agreement_K,
  width = 8,
  height = 8,
  dpi = 300
)
