# Tis script can be run after the scripts synthetic_size.R and synthetic_bands.R have been run. 
# It will join the two resulting plots in a single image.

library(patchwork)

# Output from synthetic_bands.R
p_runtime

# Output from synthetic_size.R
p_scale_runtime

p_runtime + p_scale_runtime
 
