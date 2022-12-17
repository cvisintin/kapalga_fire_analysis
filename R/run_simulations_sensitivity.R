# Load required packages
library(steps)
library(raster)
library(viridis)
library(future)
library(foreach)

# Source functions
source("R/custom_functions.R")

# Load parameter table
params <- read.csv("data/parameters_sensitivity.csv", stringsAsFactors = FALSE)

# Load spatial data
load(file = "data/spatial_rasters_background")
load(file = "data/spatial_rasters_fire_sizes")

# Create function to prepare inputs and run simulation
run_sim <- function(habitat_suitability,
                    fire_pattern,
                    juvenile_density,
                    adult_density,
                    juvenile_survival,
                    adult_survival,
                    recruitment,
                    transition,
                    juvenile_survival_sd,
                    adult_survival_sd,
                    recruitment_sd,
                    transition_sd,
                    max_ind,
                    juvenile_dispersal_dist,
                    adult_dispersal_dist,
                    juvenile_dispersal_prop,
                    adult_dispersal_prop,
                    early_surv,
                    late_surv,
                    early_trans,
                    late_trans,
                    early_recr,
                    late_recr,
                    clusters,
                    reps,
                    quoll = FALSE
) {
  
  # Create transition matrices
  trans_mat <- matrix(c(juvenile_survival, recruitment, transition, adult_survival),
                      nrow = 2, ncol = 2, byrow = TRUE)
  colnames(trans_mat) <- rownames(trans_mat) <- c('Juvenile', 'Adult')
  
  # Stochasticity matrix
  stoch_mat <- matrix(c(juvenile_survival_sd, recruitment_sd, transition_sd, adult_survival_sd),
                      nrow = 2, ncol = 2, byrow = TRUE)
  colnames(stoch_mat) <- rownames(stoch_mat) <- colnames(trans_mat)
  
  Rmax_pop <- abs(eigen(trans_mat)$values[1])
  
  # Create initial populations
  total_juveniles <- round(juvenile_density * sum(habitat_suitability[!is.na(habitat_suitability)]), 0)
  total_adults <- round(adult_density * sum(habitat_suitability[!is.na(habitat_suitability)]), 0)
  
  init_densities <- stack(habitat_suitability, habitat_suitability)
  names(init_densities) <- colnames(trans_mat)
  
  init_densities[!is.na(init_densities)] <- 0
  
  init_densities[[1]][sampleRandom(init_densities[[1]], total_juveniles, cells = TRUE)[, 1]] <- 1
  init_densities[[2]][sampleRandom(init_densities[[2]], total_adults, cells = TRUE)[, 1]] <- 1
  
  # Create carrying capacity
  k <- habitat_suitability
  k[!is.na(k)] <- max_ind
  
  ### Setup STEPS simulation ###
  
  # Specify landscape object
  landscape_obj <- landscape(population = init_densities,
                             suitability = habitat_suitability,
                             carrying_capacity = k,
                             "early_fires" = get(paste0(fire_pattern, "_early_fires")),
                             "late_fires" = get(paste0(fire_pattern, "_late_fires")))
  
  # Specify population dynamics
  pop_dyn_obj <- population_dynamics(change = growth(trans_mat,
                                                     transition_function = list(competition_density(),
                                                                                modified_transition_custom(fire_stack_early = "early_fires",
                                                                                                           fire_stack_late = "late_fires",
                                                                                                           early_surv_mult = early_surv,
                                                                                                           late_surv_mult = late_surv,
                                                                                                           early_tran_mult = early_trans,
                                                                                                           late_tran_mult = late_trans,
                                                                                                           early_recr_mult = early_recr,
                                                                                                           late_recr_mult = late_recr)),
                                                     global_stochasticity = stoch_mat),
                                     dispersal = cellular_automata_dispersal(max_cells = c(juvenile_dispersal_dist, adult_dispersal_dist), 
                                                                             use_suitability = FALSE,
                                                                             dispersal_proportion = density_dependence_dispersing(
                                                                               maximum_proportions = c(juvenile_dispersal_prop,
                                                                                                       adult_dispersal_prop))),
                                     density_dependence = ceiling_density(stages = c(1, 2)))
  
  if(quoll == TRUE) {
    pop_dyn_obj <- population_dynamics(change = growth(trans_mat,
                                                       transition_function = list(competition_density(),
                                                                                  modified_transition_custom(fire_stack_early = "early_fires",
                                                                                                             fire_stack_late = "late_fires",
                                                                                                             early_surv_mult = early_surv,
                                                                                                             late_surv_mult = late_surv,
                                                                                                             early_tran_mult = early_trans,
                                                                                                             late_tran_mult = late_trans,
                                                                                                             early_recr_mult = early_recr,
                                                                                                             late_recr_mult = late_recr,
                                                                                                             quoll = TRUE)),
                                                       global_stochasticity = stoch_mat),
                                       dispersal = cellular_automata_dispersal(max_cells = c(juvenile_dispersal_dist, adult_dispersal_dist), 
                                                                               use_suitability = FALSE,
                                                                               dispersal_proportion = density_dependence_dispersing(
                                                                                 maximum_proportions = c(juvenile_dispersal_prop,
                                                                                                         adult_dispersal_prop))),
                                       density_dependence = ceiling_density(stages = c(1, 2)))
  }
  
  # Run STEPS simulation
  
  plan(multisession, workers = reps)
  
  pop_data_totals <- foreach (i = seq_len(clusters), .combine = rbind) %do% {
    print(i)
    sim_obj <- simulation(landscape = landscape_obj,
                          population_dynamics = pop_dyn_obj,
                          demo_stochasticity = "none",
                          timesteps = 126,
                          replicates = reps,
                          verbose = FALSE,
                          future.globals = list(modified_transition_custom = modified_transition_custom,
                                                trans_mat = trans_mat,
                                                stoch_mat = stoch_mat,
                                                early_surv = early_surv,
                                                late_surv = late_surv,
                                                early_trans = early_trans,
                                                late_trans = late_trans,
                                                early_recr = early_recr,
                                                late_recr = late_recr,
                                                juvenile_dispersal_dist = juvenile_dispersal_dist,
                                                adult_dispersal_dist = adult_dispersal_dist,
                                                juvenile_dispersal_prop = juvenile_dispersal_prop,
                                                adult_dispersal_prop = adult_dispersal_prop
                          ))
    
    pop_data <- return_pop_data(sim_obj)
    rm(sim_obj)
    gc()
    pop_data
  }
  
  pop_data_totals
  
}

sp_indices <- which(params$Species != "Northern_quoll")
quoll_indices <- which(params$Species == "Northern_quoll")


reps <- 50
clusters <- 2

# Run simulations for all species except quolls
sim_pops_no_quolls <- foreach(j = sp_indices) %do% {
  run_sim(habitat_suitability = habitat_suitability,
          fire_pattern = params[j , 3],
          juvenile_density = params[j , 4],
          adult_density = params[j , 5],
          juvenile_survival = params[j , 6],
          adult_survival = params[j , 7],
          recruitment = params[j , 8],
          transition = params[j , 9],
          juvenile_survival_sd = params[j , 10],
          adult_survival_sd = params[j , 11],
          recruitment_sd = params[j , 12],
          transition_sd = params[j , 13],
          max_ind = params[j , 14],
          juvenile_dispersal_dist = params[j , 15],
          adult_dispersal_dist = params[j , 16],
          juvenile_dispersal_prop = params[j , 17],
          adult_dispersal_prop = params[j , 18],
          early_surv = params[j , 19],
          late_surv = params[j , 20],
          early_trans = params[j , 21],
          late_trans = params[j , 22],
          early_recr = params[j , 23],
          late_recr = params[j , 24],
          clusters = clusters,
          reps = reps,
          quoll = FALSE)
}

# Write out raw data
save(sim_pops_no_quolls, file = "output/all_pops_all_scenarios_sensitivity_no_quolls")

# Extract final populations for all replicates and scenarios
final_pops_no_quolls <- lapply(sp_indices,
                               function(x) data.frame(rep(paste(params[x , 1], params[x , 2], "disp", params[x , 3], "fires", sep = "_"),
                                                          nrow(sim_pops_no_quolls[[1]])),
                                                      seq_len(nrow(sim_pops_no_quolls[[1]])),
                                                      sim_pops_no_quolls[[x]][ , 126]))

final_pops_no_quolls <- do.call(rbind, final_pops_no_quolls)
colnames(final_pops_no_quolls) <- c("Scenario", "Replicate", "Final_Abundance")

# Run simulations for quolls
sim_pops_quolls <- foreach(j = quoll_indices) %do% {
  run_sim(habitat_suitability = habitat_suitability,
          fire_pattern = params[j , 3],
          juvenile_density = params[j , 4],
          adult_density = params[j , 5],
          juvenile_survival = params[j , 6],
          adult_survival = params[j , 7],
          recruitment = params[j , 8],
          transition = params[j , 9],
          juvenile_survival_sd = params[j , 10],
          adult_survival_sd = params[j , 11],
          recruitment_sd = params[j , 12],
          transition_sd = params[j , 13],
          max_ind = params[j , 14],
          juvenile_dispersal_dist = params[j , 15],
          adult_dispersal_dist = params[j , 16],
          juvenile_dispersal_prop = params[j , 17],
          adult_dispersal_prop = params[j , 18],
          early_surv = params[j , 19],
          late_surv = params[j , 20],
          early_trans = params[j , 21],
          late_trans = params[j , 22],
          early_recr = params[j , 23],
          late_recr = params[j , 24],
          clusters = clusters,
          reps = reps,
          quoll = TRUE)
}

# Write out raw data
save(sim_pops_quolls, file = "output/all_pops_all_scenarios_sensitivity_quolls")

# Extract final populations for all replicates and scenarios
final_pops_quolls <- lapply(seq_len(length(quoll_indices)),
                               function(x) data.frame(rep(paste(params[quoll_indices[x] , 1], params[quoll_indices[x] , 2], params[quoll_indices[x] , 3], "fires", sep = "_"),
                                                          nrow(sim_pops_quolls[[1]])),
                                                      seq_len(nrow(sim_pops_quolls[[1]])),
                                                      sim_pops_quolls[[x]][ , 126]))

final_pops_quolls <- do.call(rbind, final_pops_quolls)
colnames(final_pops_quolls) <- c("Scenario", "Replicate", "Final_Abundance")


# Write out final abundances
write.csv(rbind(final_pops_no_quolls, final_pops_quolls), file = "output/final_pop_size_table_sensitivity.csv", row.names = FALSE)

