# Load packages
library(rgeos)
library(ggplot2)
library(dplyr)
library(raster)
library(maptools)
library(foreach)
library(grid)
library(gridExtra)
library(gridBase)
library(gridGraphics)
library(htmlTable)
library(tableHTML)
library(rasterVis)
library(viridis)

# Create folder for figures
if(!dir.exists("figs")){
  dir.create("figs")
}

# Source functions:
source("R/custom_functions.R")

# Load parameter table
params <- read.csv("data/parameters.csv", stringsAsFactors = FALSE)

# Load in spatial layers:
load(file = "data/spatial_rasters_background")
load(file = "data/spatial_rasters_fire_sizes")

# Create EMA plot & trend plots
load(file = "output/all_pops_all_scenarios_no_quolls")
load(file = "output/all_pops_all_scenarios_quolls")

all_sim_pops <- c(sim_pops_no_quolls, sim_pops_quolls)
pops_idx <- expand.grid(1:9, 1:4)

pop_data_list <- rep(list(vector(mode = "list", length = 9)), 4)
names(pop_data_list) <- params[seq(1, length(all_sim_pops), 9) , 1]
for (i in seq_len(length(all_sim_pops) / 9)) {
  for (j in 1:9) {
    pop_data_list[[i]][[j]] <- all_sim_pops[[which(pops_idx[ , 2] == i & pops_idx[ , 1] == j)]]
  }
  names(pop_data_list[[i]]) <- paste(params[1:9 , 2], "Dispersal -", params[1:9 , 3], "Fires", sep = " ")
}


png('figs/ema_comparison.png',
    pointsize = 2,
    res = 300,
    width = 3600,
    height = 2400)

plot_ema(pop_data_list,
         all_points = FALSE,
         p = 0.01)

dev.off()

ema_data <- get_ema(pop_data_list)
write.csv(ema_data, file = "ema_data.csv", row.names = FALSE)

png('figs/possum_pop_trends.png',
    pointsize = 8,
    res = 600,
    width = 3200,
    height = 3200)

pop_data_species <- pop_data_list[["Common_brushtail_possum"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(3, 3))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()


png('figs/melomys_pop_trends.png',
    pointsize = 8,
    res = 600,
    width = 3200,
    height = 3200)

pop_data_species <- pop_data_list[["Grassland_melomys"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(3, 3))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()


png('figs/bandicoot_pop_trends.png',
    pointsize = 8,
    res = 600,
    width = 3200,
    height = 3200)

pop_data_species <- pop_data_list[["Northern_brown_bandicoot"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(3, 3))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()


png('figs/quoll_pop_trends.png',
    pointsize = 8,
    res = 600,
    width = 3200,
    height = 3200)

pop_data_species <- pop_data_list[["Northern_quoll"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(3, 3))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()



####### New table to compare percentage changes in final mean estimates of population

perc_change <- foreach(i = 1:length(pop_data_list), .combine = cbind) %do% {
  
  mean_last_pops <- unlist(lapply(pop_data_list[[i]], function(x) mean(apply(x, 1, function(y) tail(y, 1)))))
  outer(1:4, 1:4, FUN = function(x, y) (mean_last_pops[x] - mean_last_pops[y]) / mean_last_pops[x] * 100)
}

htmlTable(
  x = round(perc_change, 2),
  caption = paste("Table X. Percentage change of mean final populations",
                  "between scenarios for each species. Table is read as",
                  "change from column to row. Positive values indicate",
                  "increases and negative values indicate decreases."),
  cgroup = names(pop_data_list),
  n.cgroup=c(4, 4, 4),
  rnames = 1:4)

tableHTML(round(perc_change, 2),
          second_header = list(c(1, 4, 4, 4), c("", names(pop_data_list))),
)