##############
# Compute hot spot predictions from MAVG model
# Nil-Jana Akpinar, 2020
##############

library(spatstat)
library(pryr)
library(cubature)
library(MASS)
library(ggplot2)
library(scales)
library(assertthat)
library(svMisc)
library(tidyr)
library(future.apply)
library(VGAM)
library(grid)
library(gridExtra)

set.seed(0)

source('utils/real_bogota_utils.R')

t0 = Sys.time()


###############
# Adjustable parameters
###############

version = 1
output_folder = paste0('version_', version)

# Time steps we want predictions for
pred_ts = 2001:2190

# Parameters from data -- 00_simulate_data.R
time_steps = 6 * 365
assert_that(time_steps/(365/2) == round(time_steps/(365/2)))
cell_length = 1


###############
# Read and adjust meta data
###############

# Read meta info: district, population, victimization (per semester), report rate (per semester)
meta_path = 'metadata/bogota_victimization.csv'
meta = read.csv(meta_path)

# Scale population down
meta$population_scaled = meta$Population / 40

# Number of crimes per semester (assume victimization equals exactly one crime)
meta$n_crimes = round(meta$population_scaled * meta$Victimization)

# Add thinning rate
meta$thinning_rate = 1 - meta$Percent_reported


###############
# Start iteration
###############

dir.create(paste0('output/', output_folder, '/mavg_predictions/'), showWarnings = FALSE)

for (thinning in c(TRUE, FALSE)){
  
  cat('Start thinning =', thinning, '...\n')
  
  # Create folder for results
  dir.create(paste0('output/', output_folder, '/mavg_predictions/thinning_', thinning, '/'), showWarnings = FALSE)
  
  # Read integration data from SEPP
  if (thinning == TRUE){
    integration_data_path = paste0('output/', output_folder, '/integration_data/integration_thinned_victim_data_v', version, '_ts', time_steps, '.csv')
  } else {
    integration_data_path = paste0('output/', output_folder, '/integration_data/integration_full_victim_data_v', version, '_ts', time_steps, '.csv')
  }
  
  data = read.csv(integration_data_path)
  
  # Read parameter fit
  parameter_path = paste0('output/', output_folder, '/mavg_parameters/learned_parameter_thinning_', thinning, '.csv')
  parameters = read.csv(parameter_path)
  estim_param = parameters[1,1]
  
  
  ###############
  # Get predictions
  ###############
  
  # Number of crimes per day per cell (and add zero counts)
  bogota_mask = get_bogota_mask(cell_length = cell_length)
  bogota_mask_help = bogota_mask[bogota_mask$district != 0,]
  data = add_mask_cells(data, cell_length = cell_length)
  data$count = 1
  crime_counts = aggregate(list(count = data$count), by = list(t = data$t, cell_id = data$cell_id), FUN = sum)
  all_combis = expand.grid(t = min(data$t):max(data$t), cell_id = bogota_mask_help$cell_id)
  crime_counts = merge(x = crime_counts, y = all_combis, by = c('cell_id', 't'), all.y = TRUE)
  crime_counts$count[is.na(crime_counts$count)] = 0
  
  # Sort by cell_id and then by t
  crime_counts = crime_counts[with(crime_counts, order(cell_id, t)), ]
  
  # Add cell_wise predictions with given alpha
  add_cell_wise_pred = function(x){
    x %>% mutate(pred = as.vector(ses(count, alpha = estim_param, h = 1, initial = 'simple')$fitted))
  }
  
  crime_counts = as.data.frame(crime_counts %>%
                                 group_split(cell_id) %>%
                                 purrr::map_dfr(add_cell_wise_pred))
  
  
  # Save to file
  for (ts in pred_ts){
    rel_counts = crime_counts[crime_counts$t == ts,]
    rel_counts = merge(x = bogota_mask_help, y = rel_counts, by = c('cell_id'), all.x = TRUE)
    
    pred_path = paste0('output/', output_folder, '/mavg_predictions/thinning_', thinning, '/grid_pred_ts', ts, '.csv')
    write.csv(rel_counts, file = pred_path, row.names = FALSE)
  }
  
}
    
    t1 = Sys.time()
    time_delta = t1 - t0
    cat(paste0('Done. took: ', time_delta, ' ', units(time_delta), '\n'))
    
  