##################
# Compute integrated intensities of fitted model for evaluation period
# Nil-Jana Akpinar, 2020
##################

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

set.seed(0)

source('utils/real_bogota_utils.R')


###############
# Adjustable parameters
###############

version = 1
output_folder = paste0('version_', version)

# Time steps we want to integrate
integrate_ts = 2000:2190

# Tol for integrals
tol = 1e-5

# Parameters from data simulation
time_steps = 365 * 6
assert_that(time_steps/(365/2) == round(time_steps/(365/2)))
background_sd = 15

# Victimization data
victim_data_path = paste0('data/real_bogota_victim_data_v', version, '_ts', time_steps, '.csv')

# Time steps for training
ignore_first_t = 500 
training_time_cutoff = 4 * 500
cell_length = 1

for (ts in integrate_ts){
  if ((ts <= ignore_first_t) | (ts >= (time_steps + 1))){
    stop('Time step ts = ', ts, ' too small or large to integrate with given data.')
  }
}


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
# Compute intensities
###############

# Load full data for time steps beyond training -- no thinning
full_victim_data = read.csv(victim_data_path)
full_victim_data = full_victim_data[full_victim_data$t > ignore_first_t,]

# Get thinned data
training_data_thinned = read.csv(paste0('output/', output_folder, '/training_data/training_data_thinning_TRUE.csv'))
maxt = max(training_data_thinned$t)
candidates = full_victim_data[full_victim_data$t > maxt, c('id', 'x', 'y', 't', 'parent_index', 'district')]
thin_data = data.frame(matrix(ncol = 6, nrow = 0))

for (distr in as.character(meta$District)) {
  reportance_rate = meta[meta$District == distr, 'Percent_reported']
  rel_data = candidates[candidates$district == distr,]
  rownames(rel_data) = NULL
  bin = rbinom(n = nrow(rel_data), size = 1, prob = reportance_rate)
  thin_data = rbind(thin_data, rel_data[which(bin == 1),])
}

thinned_victim_data = rbind(training_data_thinned[,c('id', 'x', 'y', 't', 'parent_index', 'district')], thin_data)
rownames(thinned_victim_data) = NULL

# Save data
dir.create(paste0('output/', output_folder, '/integration_data/'), showWarnings = FALSE)
write.csv(thinned_victim_data, paste0('output/', output_folder, '/integration_data/integration_thinned_victim_data_v', version, '_ts', time_steps, '.csv'), row.names = FALSE)
write.csv(full_victim_data, paste0('output/', output_folder, '/integration_data/integration_full_victim_data_v', version, '_ts', time_steps, '.csv'), row.names = FALSE)

# Create folder for computed integrals
main_integral_path = paste0('output/', output_folder, '/estim_grid_integrals/')
dir.create(main_integral_path, showWarnings = FALSE)

for (thinning in c(TRUE, FALSE)){
  
  cat('Start thinning =', thinning, '...\n')
  
  # Read parameter fits
  parameter_path = paste0('output/', output_folder, '/learned_parameters/learned_parameters_thinning_', thinning, '.csv')
  parameters = read.csv(parameter_path)
  n_iter = nrow(parameters)
  estim_mu_bar = parameters$mu_bar[n_iter]
  estim_theta = parameters$theta[n_iter]
  estim_omega = parameters$omega[n_iter]
  estim_sigma_x = parameters$sigma_x[n_iter]
  estim_sigma_y = parameters$sigma_y[n_iter]
  
  # Get cells to integrate
  mask = get_bogota_mask(cell_length = cell_length)
  bogota_mask_help = mask[mask$district != 0,]
  
  # Get reasonable truncation parameter
  f = function(N, data, theta, omega, sigma_x, sigma_y){
    nrow(data) * theta * omega * exp(-omega * N) * 2 * pi * sigma_x * sigma_y
  }
  
  if (thinning == TRUE) {
    
    full_int_600 = read.csv('output/version_1/estim_grid_integrals/thinning_TRUE/full_integrals/grid_integrals_ts2000.csv') # Precomputed example version
    sp = sort(full_int_600$int, TRUE)[ceiling(nrow(bogota_mask_help) * 0.2)]
    max_absolute_error = sp/100
    fN = function(N) pryr::partial(f, data = thinned_victim_data, theta = estim_theta, omega = estim_omega, sigma_x = estim_sigma_x, sigma_y = estim_sigma_y)(N)
    Ns = 1:1000
    fs = unlist(lapply(Ns, fN))
    N = which(fs <= max_absolute_error)[1]
    cat(paste0('Integrate intensities considering crimes of the last N = ', N, ' days.\n'))
  
    } else if (thinning == FALSE) {
    
      full_int_600 = read.csv('output/version_1/estim_grid_integrals/thinning_FALSE/full_integrals/grid_integrals_ts2000.csv') # Precomputed example version
      sp = sort(full_int_600$int, TRUE)[ceiling(nrow(bogota_mask_help) * 0.2)]
      max_absolute_error = sp/100
      fN = function(N) pryr::partial(f, data = full_victim_data, theta = estim_theta, omega = estim_omega, sigma_x = estim_sigma_x, sigma_y = estim_sigma_y)(N)
      Ns = 1:1000
      fs = unlist(lapply(Ns, fN))
      N = which(fs <= max_absolute_error)[1]
      cat(paste0('Integrate intensities considering crimes of the last N = ', N, ' days.\n'))
      
  } else {
    stop('No valid value for thinning.')
  }

  # Integrate background intensity
  res = get_background_integrals(centers = data.frame(x = 0, y = 0),
                                 cell_length = cell_length,
                                 total_mu_bar = estim_mu_bar,
                                 background_sd = background_sd)
  
  colnames(res) = c('x', 'y', 'district', 'cell_id', 'background_int')
  bogota_mask_help = merge(x = bogota_mask_help, y = res, all.x = TRUE)
  
  # Get fit of aftershock intensity function
  estim_aftershock_intensity = function(x, y, t) pryr::partial(aftershock_intensity,
                                                               x_center = 0,
                                                               y_center = 0,
                                                               t_center = 0,
                                                               sigma_x = estim_sigma_y,
                                                               sigma_y = estim_sigma_y,
                                                               omega = estim_omega,
                                                               theta = estim_theta)(x,y,t)
  # Integrate aftershock intensities
  dir.create(paste0(main_integral_path,'thinning_', thinning, '/'), showWarnings = FALSE)
  
  # Compute aftershock intensities
  plan(multiprocess)
  YY = future_lapply(integrate_ts, function(ts){
    
    cat(paste0('Compute integrals for t = ', ts, '...'))
    t0 = Sys.time()
    
    # Path to save computed integrals
    integrals_path = paste0(main_integral_path,'thinning_', thinning, '/grid_integrals_ts', ts, '.csv')
    
    # Get relevant data to this time step
    if (thinning == TRUE) {
      if ((ts <= min(thinned_victim_data$t)) | (ts > (max(thinned_victim_data$t) + 1))){
        stop(paste0('Time step ts = ', ts, ' is too large or small to integrate with the given data.'))
      } else {
        rel_data = thinned_victim_data[(thinned_victim_data$t < ts) & (ts - thinned_victim_data$t < N),]
      }
    } else if (thinning == FALSE) {
      if ((ts <= min(full_victim_data$t)) | (ts > (max(full_victim_data$t) + 1))){
        stop(paste0('Time step ts = ', ts, ' is too large or small to integrate with the given data.'))
      } else {
        rel_data = full_victim_data[(full_victim_data$t < ts) & (ts - full_victim_data$t < N),]
      }
    }
    
    res = get_aftershock_integrals(centers = rel_data,
                                   ts = ts,
                                   sds = c(estim_sigma_x, estim_sigma_y),
                                   theta = estim_theta, 
                                   omega = estim_omega, 
                                   cell_length = cell_length)
    
    colnames(res) = c('x', 'y', 'district', 'cell_id', 'aftershock_int')
    
    res = merge(x = bogota_mask_help, y = res, all.x = TRUE)
    res$int = res$background_int + res$aftershock_int
    
    # Save to file
    write.csv(res, file = integrals_path, row.names = FALSE)
    
    t1 = Sys.time()
    time_delta = t1 - t0
    cat(paste0('Done. took: ', time_delta, ' ', units(time_delta), '\n'))
    
  })
  
}
