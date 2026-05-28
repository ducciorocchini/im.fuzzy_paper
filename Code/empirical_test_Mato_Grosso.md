# Empirical fuzzy classification on Mato Grosso imagery (imageRy) 

```r
# with 2 classes and class-wise histograms.

library(terra)
library(imageRy)
library(ggplot2)
library(ggridges)

# import
mg <- rast("~/Desktop/STScI-01EVSV9JYWYARC4A8336KMKNTQ.jpeg")   
mg <- flip(mg)

# Ensure we have a multi-band SpatRaster (JPEG typically yields 3 bands: R,G,B)
mg
nlyr(mg)

# 2) Run fuzzy classification with 4 classes
# (Assumes you have im.fuzzy() available in your environment/package)
res <- im.fuzzy(mg, num_clusters = 4, seed = 42)

mem <- res$memberships   # SpatRaster with 2 layers:
# membership to class 1 and class 2
dist <- res$distances    # SpatRaster with 2 layers: 
# distances to centroid 1 and 2

# 3) Output plots

# A) Brightness from RGB
brightness <- (mg[[1]] + mg[[2]] + mg[[3]]) / 3
names(brightness) <- "brightness"

# B) Sample brightness + memberships
mem <- res$memberships
K <- nlyr(mem)
names(mem) <- paste0("Class ", 1:K)

set.seed(42)
samp <- spatSample(
  c(brightness, mem),
  size = 200000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

x <- samp$brightness
M <- as.matrix(samp[, -1])

# C) Long format (fuzzy weights)
long <- data.frame(
  brightness = rep(x, times = K),
  class = factor(rep(colnames(M), each = length(x))),
  weight = as.vector(M)
)

# D)
ggplot(long, aes(x = brightness, weight = weight, fill = class)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  theme_minimal() +
  labs(
    x = "RGB brightness",
    y = "Weighted density",
    title = "Fuzzy brightness densities"
  )
```

