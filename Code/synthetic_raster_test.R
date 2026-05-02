library(terra)
library(cluster)      # silhouette
# install.packages(c("e1071", "mclust", "aricode", "cluster"))
# install.packages("mclust") # ARI
library(mclust)      # ARI
# install.packages("aricode") # NMI
library(aricode)     # NMI
library(e1071)
library(imageRy) # im.fuzzy function

set.seed(42)


# -----------------------
# 1. Generate synthetic raster
# -----------------------
n <- 100 * 100
K <- 4
B <- 3

# True class labels
true_labels <- sample(1:K, n, replace = TRUE)

# Generate Gaussian clusters
centers <- matrix(runif(K * B, 0, 255), ncol = B)

X <- matrix(NA, n, B)
for (i in 1:n) {
  X[i, ] <- rnorm(B, mean = centers[true_labels[i], ], sd = 20)
}

# Convert to raster
r <- rast(nrows = 100, ncols = 100, nlyrs = B)
values(r) <- X

# -----------------------
# 2. Methods
# -----------------------

# k-means
km <- kmeans(X, centers = K)
km_labels <- km$cluster

# fuzzy c-means
fcm <- cmeans(X, centers = K, m = 2)
fcm_labels <- apply(fcm$membership, 1, which.max)

# im.fuzzy (your function)
res <- im.fuzzy(r, num_clusters = K, seed = 42, do_plot = FALSE)
mem <- values(res$memberships)
imf_labels <- apply(mem, 1, which.max)

# -----------------------
# 3. External validation
# -----------------------

acc <- function(true, pred) mean(true == pred)

results <- data.frame(
  Method = c("k-means", "fuzzy c-means", "im.fuzzy"),
  Accuracy = c(
    acc(true_labels, km_labels),
    acc(true_labels, fcm_labels),
    acc(true_labels, imf_labels)
  ),
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
