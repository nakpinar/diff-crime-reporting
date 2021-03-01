##############
# Integrate true intensities from data (not dependent on fitted model; only recompute when data is updated)
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

source('utils/real_bogota_utils.R')
source('utils/model_utils.R')

set.seed(0)


###############
# Adjustable parameters
###############

# Time steps we want to integrate; only necessary for evaluation period
integrate_ts = 500:2190

# Tolerance for integrals
tol = 1e-5

# Parameters for true superset data -- same as in 00_simulate_data.R
version = 1
cell_length = 1
time_steps = 365 * 6
assert_that(time_steps/(365/2) == round(time_steps/(365/2)))

mu_bar = 100 
theta = 15
omega = 0.2
sigma_x = 0.1
sigma_y =  0.1 

data_name = paste0('real_bogota_data_v', version, '_ts', time_steps)
superset_data_file = paste0('data/real_bogota_superset_data_v', version, '_ts', time_steps, '.csv')
superset_data = read.csv(superset_data_file)


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

# Add thinning rates for each districts that give us the desired victimization in each district -- based on superset data
n_semesters = time_steps/(365/2)
meta$keep_rate = 0
cat('Compute keep rates according to victimization...\n')
for (distr in as.character(meta$District)){
  desired_n_crimes = meta[meta$District == distr, 'n_crimes'] * n_semesters
  sampled_n_crimes = nrow(superset_data[superset_data$district == distr,])
  if (desired_n_crimes > sampled_n_crimes){
    stop(paste0('Not enough crimes for ', distr, '.\nSupersample: ', sampled_n_crimes, '\nDesired: ', desired_n_crimes))
  } else {
    cat('For ', distr, ': sampled_n = ', sampled_n_crimes, ', desired_n = ', desired_n_crimes, '\n')
  }
  meta[meta$District == distr, 'keep_rate'] = desired_n_crimes/sampled_n_crimes
}


###############
# Compute intensities
###############

# Start integration
mask = get_bogota_mask(cell_length = cell_length)
bogota_mask_help = mask[mask$district != 0,]

## Speed up computation by N truncation -- only consider events that are not too far away in time for triggering function
f = function(N){
  nrow(superset_data) * theta * omega * exp(-omega * N) * 2 * pi * sigma_x * sigma_y
}

# Smallest value we care to estimate precisely (within 1% relative error) assuming we would at most predict 20% of the city as hot spots
full_int_1000 = read.csv(paste0('data/real_bogota_data_v1_ts2190/grid_integrals_full/grid_integrals_ts1000.csv'))
sp = sort(full_int_1000$superset_int, TRUE)[ceiling(nrow(bogota_mask_help) * 0.2)]
max_absolute_error = sp/100
sp
max_absolute_error

# Find smallest N we can choose in order to make smaller errors than max_absolute error
Ns = 1:200
fs = unlist(lapply(Ns, f))
N = which(fs <= max_absolute_error)[1]
cat(paste0('Integrate intensities considering crimes of the last N = ', N, ' days.\n'))

# Create folders for intensities
dir.create(paste0('data/', data_name, '/'), showWarnings = FALSE)
dir.create(paste0('data/', data_name, '/grid_integrals/'), showWarnings = FALSE)

# Integrate background (independent of t and data)
cat('Compute true superset integral background...')
t0 = Sys.time()

centers = data.frame(do.call(rbind, get_centers('realistic_bogota')))
colnames(centers) = c('x', 'y')

res = get_background_integrals(centers = centers,
                              cell_length = cell_length,
                              total_mu_bar = mu_bar,
                              background_sd = 4.5)

bogota_mask_help = merge(x = bogota_mask_help, y = res, all.x = TRUE)

t1 = Sys.time()
time_delta = t1 - t0
cat(paste0('Done. took: ', time_delta, ' ', units(time_delta), '\n'))

# Integrate aftershock intensities
plan(multiprocess)

#for (ts in integrate_ts){
YY = future_lapply(integrate_ts, function(ts){
  
  cat(paste0('Compute integrals for t = ', ts, ':\n'))
  
  # Get paths to save computed integrals
  integrals_path = paste0('data/', data_name, '/grid_integrals/grid_integrals_ts', ts, '.csv')
  
  ########
  # Compute true superset integrals
  cat('Compute true superset integral...')
  t0 = Sys.time()
  
  # Get aftershock integrals
  rel_superset_data = superset_data[(superset_data$t < ts) & (ts - superset_data$t < N),]
  res = get_aftershock_integrals(centers = rel_superset_data,
                                ts = ts,
                                sds = c(sigma_x, sigma_y),
                                theta = theta, 
                                omega = omega, 
                                cell_length = cell_length)
  
  res = merge(x = bogota_mask_help, y = res, all.x = TRUE)
  res$superset_int = res$superset_background_int + res$superset_aftershock_int
  
  t1 = Sys.time()
  time_delta = t1 - t0
  cat(paste0('Done. took: ', time_delta, ' ', units(time_delta), '\n'))
  
  ########
  # Compute true victimization integrals from true superset integrals (add to same df)
  cat('Add victimization integrals...')
  
  res$victim_int = 0
  for (d in as.character(meta$District)){
    keep_rate = meta[meta$District == d, 'keep_rate']
    res[res$district == d, 'victim_int'] = res[res$district == d, 'superset_int'] * keep_rate
  }
  
  cat('Done.\n')
  
  # Save results
  write.csv(res, file = integrals_path, row.names = FALSE)
  
})

