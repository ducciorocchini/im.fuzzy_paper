# ============================================================
# BASELINE COMPUTATIONAL COMPARISON
#
# Comparison among:
#   1. k-means
#   2. fuzzy c-means
#   3. im.fuzzy()
#
# Metrics:
#   - Runtime (seconds)
#   - Peak RAM used during execution (MiB)
#
# Experimental design:
#   - Vary number of pixels N
#   - Keep number of bands B fixed
#   - Keep number of clusters K fixed
#   - Repeat each configuration several times
#
# IMPORTANT:
# For fairness, preprocessing is included in the timed/memory
# benchmark for all three methods.
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

library(terra)
library(e1071)
library(ggplot2)
library(imageRy)
library(peakRAM)
library(dplyr)

# If needed:
# install.packages(c("peakRAM", "dplyr"))


# ------------------------------------------------------------
# 2. Experimental parameters
# ------------------------------------------------------------

# Number of pixels to test
N_values <- c(
  1e4,
  5e4,
  1e5,
  2.5e5,
  5e5
)

# Fixed number of bands
B <- 4

# Fixed number of clusters
K <- 4

# Fuzzifier
m <- 2

# Number of repetitions
n_rep <- 10

# Base random seed
base_seed <- 42


# ------------------------------------------------------------
# 3. Storage object
# ------------------------------------------------------------

results <- data.frame(
  N = integer(),
  Replicate = integer(),
  Method = character(),
  Runtime = numeric(),
  Peak_RAM_MiB = numeric(),
  stringsAsFactors = FALSE
)


# ============================================================
# 4. Helper function
# Scale each band independently to 0-255
# ============================================================

scale_255 <- function(x) {

  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)

  if (xmax == xmin) {
    return(rep(0, length(x)))
  }

  255 * (x - xmin) / (xmax - xmin)
}


# ============================================================
# 5. Main experiment
# ============================================================

for (N in N_values) {

  cat("\n")
  cat("============================================\n")
  cat("Requested number of pixels:", N, "\n")
  cat("============================================\n")

  # ----------------------------------------------------------
  # Raster dimensions approximately matching N
  # ----------------------------------------------------------

  nrow_r <- floor(sqrt(N))
  ncol_r <- ceiling(N / nrow_r)

  N_actual <- nrow_r * ncol_r

  cat("Actual number of pixels:", N_actual, "\n")


  # ----------------------------------------------------------
  # Repetitions
  # ----------------------------------------------------------

  for (rep in seq_len(n_rep)) {

    cat(
      "\nN =",
      N_actual,
      "- replicate",
      rep,
      "of",
      n_rep,
      "\n"
    )


    # --------------------------------------------------------
    # Generate synthetic data
    # --------------------------------------------------------

    set.seed(base_seed + rep)

    X <- matrix(
      rnorm(N_actual * B),
      nrow = N_actual,
      ncol = B
    )


    # --------------------------------------------------------
    # Build equivalent SpatRaster for im.fuzzy()
    # --------------------------------------------------------

    r <- terra::rast(
      nrows = nrow_r,
      ncols = ncol_r,
      nlyrs = B
    )

    terra::values(r) <- X


    # ========================================================
    # METHOD 1: K-MEANS
    # ========================================================

    cat("Running k-means\n")

    gc()
    set.seed(base_seed)

    ram_km <- peakRAM::peakRAM({

      # Preprocessing included in benchmark
      X_km <- apply(
        X,
        2,
        scale_255
      )

      X_km <- as.matrix(X_km)

      km <- stats::kmeans(
        X_km,
        centers = K
      )

    })


    runtime_km <- ram_km$Elapsed_Time_sec[1]

    # On this system, Peak_RAM_Used_MiB is returned in bytes.
    # Convert explicitly to MiB.
    memory_km <- ram_km$Peak_RAM_Used_MiB[1] / 1024^2


    results <- rbind(
      results,
      data.frame(
        N = N_actual,
        Replicate = rep,
        Method = "k-means",
        Runtime = runtime_km,
        Peak_RAM_MiB = memory_km
      )
    )


    rm(km, X_km)
    gc()


    # ========================================================
    # METHOD 2: FUZZY C-MEANS
    # ========================================================

    cat("Running fuzzy c-means\n")

    gc()
    set.seed(base_seed)

    ram_fcm <- peakRAM::peakRAM({

      # Preprocessing included in benchmark
      X_fcm <- apply(
        X,
        2,
        scale_255
      )

      X_fcm <- as.matrix(X_fcm)

      fcm <- e1071::cmeans(
        X_fcm,
        centers = K,
        m = m
      )

    })


    runtime_fcm <- ram_fcm$Elapsed_Time_sec[1]

    memory_fcm <- ram_fcm$Peak_RAM_Used_MiB[1] / 1024^2


    results <- rbind(
      results,
      data.frame(
        N = N_actual,
        Replicate = rep,
        Method = "fuzzy c-means",
        Runtime = runtime_fcm,
        Peak_RAM_MiB = memory_fcm
      )
    )


    rm(fcm, X_fcm)
    gc()


    # ========================================================
    # METHOD 3: im.fuzzy()
    # ========================================================

    cat("Running im.fuzzy\n")

    gc()
    set.seed(base_seed)

    ram_imf <- peakRAM::peakRAM({

      imf <- imageRy::im.fuzzy(
        r,
        num_clusters = K,
        m = m,
        seed = base_seed,
        do_plot = FALSE
      )

    })


    runtime_imf <- ram_imf$Elapsed_Time_sec[1]

    memory_imf <- ram_imf$Peak_RAM_Used_MiB[1] / 1024^2


    results <- rbind(
      results,
      data.frame(
        N = N_actual,
        Replicate = rep,
        Method = "im.fuzzy",
        Runtime = runtime_imf,
        Peak_RAM_MiB = memory_imf
      )
    )


    rm(imf)
    gc()


    # --------------------------------------------------------
    # Clean large objects before next repetition
    # --------------------------------------------------------

    rm(
      X,
      r,
      ram_km,
      ram_fcm,
      ram_imf
    )

    gc()
  }
}


# ============================================================
# 6. Inspect raw results
# ============================================================

cat("\n")
cat("============================================\n")
cat("RAW RESULTS\n")
cat("============================================\n")

print(results)

cat("\nColumn names:\n")
print(names(results))


# ============================================================
# 7. Summary statistics
# ============================================================

summary_clean <- results %>%
  group_by(N, Method) %>%
  summarise(
    Runtime_mean = mean(Runtime, na.rm = TRUE),
    Runtime_sd = sd(Runtime, na.rm = TRUE),
    Peak_RAM_mean = mean(Peak_RAM_MiB, na.rm = TRUE),
    Peak_RAM_sd = sd(Peak_RAM_MiB, na.rm = TRUE),
    .groups = "drop"
  )


cat("\n")
cat("============================================\n")
cat("SUMMARY RESULTS\n")
cat("============================================\n")

print(summary_clean)

cat("\nSummary column names:\n")
print(names(summary_clean))


# ============================================================
# 8. Runtime plot
# ============================================================

p_runtime <- ggplot(
  summary_clean,
  aes(
    x = N,
    y = Runtime_mean,
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

  geom_errorbar(
    aes(
      ymin = pmax(
        0,
        Runtime_mean - Runtime_sd
      ),
      ymax = Runtime_mean + Runtime_sd
    ),
    width = 0
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Computational runtime comparison",
    subtitle = paste(
      "B =", B,
      "bands; K =", K,
      "clusters"
    ),
    x = "Number of pixels",
    y = "Runtime (seconds)",
    color = "Method"
  )


print(p_runtime)


# ============================================================
# 9. Peak RAM plot
# ============================================================

p_memory <- ggplot(
  summary_clean,
  aes(
    x = N,
    y = Peak_RAM_mean,
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

  geom_errorbar(
    aes(
      ymin = pmax(
        0,
        Peak_RAM_mean - Peak_RAM_sd
      ),
      ymax = Peak_RAM_mean + Peak_RAM_sd
    ),
    width = 0
  ) +

  theme_minimal(
    base_size = 14
  ) +

  labs(
    title = "Peak memory comparison",
    subtitle = paste(
      "B =", B,
      "bands; K =", K,
      "clusters"
    ),
    x = "Number of pixels",
    y = "Peak RAM used (MiB)",
    color = "Method"
  )


print(p_memory)


# ============================================================
# 10. Results at largest image size
# ============================================================

largest_N <- max(
  summary_clean$N
)

largest_results <- summary_clean %>%
  filter(N == largest_N)


cat("\n")
cat("============================================\n")
cat("RESULTS AT LARGEST IMAGE SIZE\n")
cat("============================================\n")

print(largest_results)


# ============================================================
# 11. Save numerical results
# ============================================================

write.csv(
  results,
  "baseline_computational_raw.csv",
  row.names = FALSE
)

write.csv(
  summary_clean,
  "baseline_computational_summary.csv",
  row.names = FALSE
)

write.csv(
  largest_results,
  "baseline_computational_largestN.csv",
  row.names = FALSE
)


# ============================================================
# 12. Save figures
# ============================================================

ggsave(
  filename = "baseline_runtime_comparison.png",
  plot = p_runtime,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "baseline_peakRAM_comparison.png",
  plot = p_memory,
  width = 8,
  height = 5,
  dpi = 300
)

