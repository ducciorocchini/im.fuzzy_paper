# Output maps

For visual comparison, clustering outputs obtained with the three methods are provided:

1. **k-means**  
   Crisp clustering partition, in which each pixel is assigned exclusively to one cluster.

   File: k_means.png

3. **Fuzzy C-means (FCM)**  
   Fuzzy membership maps showing the degree of membership of each pixel to each cluster.

   File: fuzzy_c_means.png

5. **`im.fuzzy()`**  
   Fuzzy membership maps derived from pixel-to-centroid distances using k-means-derived cluster centers.

   File: monsummano_imfuzzy_memberships_K2.png

For the fuzzy methods, one membership map is provided for each cluster, with membership values ranging from 0 to 1. Cluster identities were matched across methods to facilitate visual comparison of corresponding spatial patterns.
