library(terra)
library(dplyr)

# 1. Load your cluster raster
cluster_raster <- rast("forestcover_grid_sf.tif")

# 2. Identify patches (clumps) of same cluster (8-direction adjacency)
patch_raster <- patches(cluster_raster, directions = 8, values = TRUE)

# 3. Bind the two rasters into one SpatRaster (so they keep their names)
combo <- c(patch_raster, cluster_raster)

# 4. Convert to a data.frame, dropping all NA cells at once
#    This will give you two columns named after the layers
df <- as.data.frame(combo, na.rm = TRUE)

# 5. Rename columns for clarity
#    (adjust these names if your layer-names differ)
names(df) <- c("patch_id", "cluster")

# 6. Now compute patch-level stats
patch_stats <- df %>%
  group_by(patch_id, cluster) %>%
  summarise(
    n_pixels = n(),
    area_ha  = n_pixels * 0.01,  # 100 m² per pixel → 0.01 ha
    .groups = "drop"
  )

# 7. Summarize per cluster
cluster_stats <- patch_stats %>%
  group_by(cluster) %>%
  summarise(
    n_patches               = n(),
    mean_area_ha            = mean(area_ha),
    sd_area_ha              = sd(area_ha),
    patches_above_1ha       = sum(area_ha > 1),
    pct_patches_above_1ha   = 100 * patches_above_1ha / n_patches,
    mean_area_ge1ha         = mean(area_ha[area_ha >= 1], na.rm = TRUE),
    total_area_ha           = sum(area_ha),
    area_above_1ha          = sum(area_ha[area_ha > 1]),
    pct_area_above_1ha      = 100 * area_above_1ha / total_area_ha,
    .groups = "drop"
  ) 

# 8. View or export
print(cluster_stats)
#write.csv(cluster_stats, "cluster_summary_window10_k10_wardD2.csv", row.names = FALSE)
#write_xlsx(cluster_stats, "cluster_summary_window10_k10_wardD2.xlsx")
