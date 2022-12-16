# Load packages
library(steps)
library(raster)
library(viridis)
library(future)
library(foreach)

# Source functions
source("R/custom_functions.R")

# Load in spatial layers
load(file = "data/spatial_rasters_background")

# Get index of non-NA cells
idx_no_na <- which(!is.na(habitat_suitability[]))

# Create table of number of cells burned and proportions early versus late
burn_table <- data.frame("timestep" = 1:126, "early_fire" = NA, "late_fire" = NA, "proportion" = NA)

# Create vectors of early and late timesteps
early_fire_timesteps <- seq(3, 126, 6)
late_fire_timesteps <- seq(5, 126, 6)

# Calculate early fire proportions
for (i in early_fire_timesteps) {
  burn_table[i, 2] <- sum(recorded_early_fires[[i]][idx_no_na] == 1)
  burn_table[i, 4] <- burn_table[i, 2] / length(idx_no_na)
}

# Calculate late fire proportions
for (i in late_fire_timesteps) {
  burn_table[i, 3] <- sum(recorded_late_fires[[i]][idx_no_na] == 1)
  burn_table[i, 4] <- burn_table[i, 3] / length(idx_no_na)
}

# Create new empty raster mask
mask <- habitat_suitability
mask[which(mask[] == 1)] <- 0


### Dispersed random 1-hectare fires ###

# Create empty stack for dispersed fires
dispersed_early_fires <- brick(replicate(126, mask))
dispersed_late_fires <- brick(replicate(126, mask))

# Sample and replace cells with fires in early timesteps
for (i in early_fire_timesteps) {
  set.seed(i)
  early_idx <- sample(idx_no_na, burn_table[i, 2])
  dispersed_early_fires[[i]][early_idx] <- 1 
}
names(dispersed_early_fires) <- paste0(1:nlayers(dispersed_early_fires))

# Sample and replace cells with fires in late timesteps that were
# not burnt in early timesteps
for (i in late_fire_timesteps) {
  early_idx <- which(dispersed_early_fires[[i - 2]][] == 1)
  late_idx <- setdiff(idx_no_na, early_idx)
  set.seed(i + 100)
  late_idx <- sample(late_idx, burn_table[i, 3])
  dispersed_late_fires[[i]][late_idx] <- 1 
}
names(dispersed_late_fires) <- paste0(1:nlayers(dispersed_late_fires))


### Clumped random large fires ###

# Create empty stack for clumped fires
clumped_early_fires <- brick(replicate(126, mask))
clumped_late_fires <- brick(replicate(126, mask))

candidate_cells <- c(head(idx_no_na, 5000), tail(idx_no_na, 5000))
candidate_cells <- sample(candidate_cells)

# Sample and replace cells with fires in early timesteps
for (i in early_fire_timesteps) {
  print(paste0("Timestep ", i))
  set.seed(i + 300)
  start_cell <- sample(candidate_cells, 1)
  print(start_cell)
  clumped_early_fires[[i]] <- create_fires(input_raster = habitat_suitability,
                                           masked_cells = NULL,
                                           start_cell = start_cell,
                                           n_cells_burned = burn_table[i, 2],
                                           cell_radius = 3) 
}
names(clumped_early_fires) <- paste0(1:nlayers(clumped_early_fires))
plot(clumped_early_fires[[early_fire_timesteps]], maxnl = 24)
plot(sum(clumped_early_fires))

# Check what percentage the variances are from the original totals
sim_early_burned_cells <- sapply(early_fire_timesteps, FUN =  function(x) sum(clumped_early_fires[[x]][], na.rm = TRUE)) 
(error_early_burned_cells <- abs((sim_early_burned_cells - burn_table[early_fire_timesteps, 2]) / burn_table[early_fire_timesteps, 2]) * 100)

# Sample and replace cells with fires in late timesteps that were not burnt in early timesteps
for (i in late_fire_timesteps) {
  print(paste0("Timestep ", i))
  early_idx <- which(clumped_early_fires[[i - 2]][] == 1)
  set.seed(i + 500)
  start_cell <- sample(setdiff(idx_no_na, early_idx), 1)
  print(start_cell)
  clumped_late_fires[[i]] <- create_fires(input_raster = habitat_suitability,
                                          masked_cells = early_idx,
                                          start_cell = start_cell,
                                          n_cells_burned = burn_table[i, 3],
                                          cell_radius = 3)
}
names(clumped_late_fires) <- paste0(1:nlayers(clumped_late_fires))
plot(clumped_late_fires[[late_fire_timesteps]], maxnl = 24)
plot(sum(clumped_late_fires))

# Check what percentage the variances are from the original totals
sim_late_burned_cells <- sapply(late_fire_timesteps, FUN =  function(x) sum(clumped_late_fires[[x]][], na.rm = TRUE)) 
(error_late_burned_cells <- abs((sim_late_burned_cells - burn_table[late_fire_timesteps, 3]) / burn_table[late_fire_timesteps, 3]) * 100)

# Write out spatial fire data to be read back into R simulations
save(dispersed_early_fires,
     dispersed_late_fires,
     clumped_early_fires,
     clumped_late_fires,
     file = "data/spatial_rasters_fire_sizes",
     compress = TRUE)
