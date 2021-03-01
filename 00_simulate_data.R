###############
# Simulate synthetic Bogota crime data according to surveyed victimization and victim crime reporting rates
# Nil-Jana Akpinar, 2020
###############

library(spatstat)
library(pryr)
library(cubature)
library(MASS)
library(ggplot2)
library(scales)
library(assertthat)
library(svMisc)
library(dplyr)
library(future.apply)

source('utils/model_utils.R')
source('utils/real_bogota_utils.R')

set.seed(0)

t0 = Sys.time()


###############
# Parameters
###############

version = 1
time_steps = 365 * 6
assert_that(time_steps/(365/2) == round(time_steps/(365/2)))

# Data parameters for high intensity SEPP
mu_bar = 100
theta = 15
omega = 0.2
sigma_x = 0.1
sigma_y =  0.1 


###############
# Read and adjust meta data
###############

# Read meta info: district, population, victimization (per semester), report rate (per semester)
meta_path = 'metadata/bogota_victimization.csv'
meta = read.csv(meta_path)

# Scale population down
meta$population_scaled = meta$Population / 40

# Number of crimes per semester
meta$n_crimes = round(meta$population_scaled * meta$Victimization)

# Add thinning rate
meta$thinning_rate = 1 - meta$Percent_reported


###############
# Superset data -- background
###############

t00 = Sys.time()

# Data paths
superset_data_file = paste0('../data/realistic_bogota/superset_data_v', version, '_ts', time_steps, '.csv')
data_file = paste0('../data/realistic_bogota/data_v', version, '_ts', time_steps, '.csv')

# Get partial background functions
centers = get_centers('realistic_bogota')
mu_partial = mu_bar/length(centers)
for (ii in 1:length(centers)){
  local({
    idx = ii
    assign(paste0('background_intensity_part', idx), get_unimodal_background_intensity_function(x, y, centers[[idx]], mu_partial), pos = .GlobalEnv)
  })
}

# Sample background events for each time step
data = data.frame(matrix(ncol = 5, nrow = 0))
xrange = c(-13.15241, 13.15241)
yrange = c(-31.35409, 31.35409)

for (t in 1:time_steps){
  for (ii in 1:length(centers)){
    pp = rpoispp(lambda = get(paste0('background_intensity_part', ii)), 1, win = owin(xrange,yrange))
    data = rbind(data, cbind(pp$x, pp$y, rep(t,pp$n), rep(0,pp$n), rep(ii, pp$n)))
  }
}

colnames(data) = c('x','y','t','parent_index','group')
data$group = as.factor(data$group)

# Add id column
rownames(data) = NULL
data$id = 1:nrow(data)

# Plot background data for first month as check
# ggplot() + 
#   geom_bogota() + 
#   theme_bogota() + 
#   geom_point(data = data[data$t <= 30,], aes(x = x, y = y, color = group))

# Remove crime outside city bounds
data = add_districts(data, cell_length = 1)
data = data[data$district != 0,]

# Plot superset background crime in first time steps
# ggplot() + 
#   geom_bogota() + 
#   theme_bogota() + 
#   add_grid() +
#   geom_point(data = data[data$t <= 100,], aes(x = x, y = y, color = district))

t11 = Sys.time()
t_del = t11 - t00
cat(paste0('Superset background done. Took ', t_del, ' ', units(t_del), '\n'))


###############
# Superset data -- aftershock
###############

t00 = Sys.time()

# Delete group column
data = data[,c('id', 'x', 'y', 't', 'parent_index', 'district')]
background_data = data

# Number of iterations and cumulative number of background events after iterations 
l = 0
n_background = c(nrow(data))

cat('Sample aftershock events.\n')

plan(multiprocess)

# Expected number of offspring events per crime
m = theta*2*pi*sigma_x*sigma_y

while (nrow(background_data) > 0){
  
  l = l + 1
  
  # Add number of offspring events
  background_data$n_offspring = rpois(nrow(background_data), m)
  sum_offspring = sum(background_data$n_offspring)
  
  cat(paste0('Iteration: ', l, '     Added aftershocks: ', sum_offspring, '\n'))
  
  # Draw uncentered offspring samples, we can draw each coordinate separately
  x = rnorm(sum_offspring, mean = 0, sd = sigma_x)
  y = rnorm(sum_offspring, mean = 0, sd = sigma_y)
  t = ceiling(rexp(sum_offspring, rate = omega))
  
  # Initialize offspring data frame
  offspring_data = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(offspring_data) = c('id', 'x', 'y', 't', 'parent_index')

  YY = future_lapply(1:nrow(background_data), function(ii){
    
    if (background_data[ii, 'n_offspring'] > 0){
      
      x_center = background_data[ii, 'x']
      y_center = background_data[ii, 'y']
      t_center = background_data[ii, 't']
      parent_id = background_data[ii, 'id']
      n_offspring = background_data[ii, 'n_offspring']
    
      # Get valid indices from x, y, t
      start_idx = cumsum(background_data$n_offspring)[ii] - n_offspring
      data_offspring_partial = data.frame(matrix(ncol = 4, nrow = 0))
      
      for (j in 1:n_offspring){
        data_offspring_partial = rbind(data_offspring_partial, c(x_center + x[start_idx + j], y_center + y[start_idx + j], t_center + t[start_idx + j], parent_id))
      }
      
      colnames(data_offspring_partial) = c('x', 'y', 't', 'parent_index')
      data_offspring_partial
      
        
    } else {
        
      data_offspring_partial = data.frame(matrix(ncol = 4, nrow = 0))
      colnames(data_offspring_partial) = c('x', 'y', 't', 'parent_index')
      data_offspring_partial
      
      }
  })
  
  # Assign unique indices
  offspring_data = bind_rows(YY, .id = NULL)
  if (nrow(offspring_data) > 0){
    offspring_data$id = (max(data$id) + 1):(max(data$id) + nrow(offspring_data))
  } else {
    offspring_data = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(offspring_data) = c('x', 'y', 't', 'parent_index', 'id')
  }
  
  # Erase out of bounds data
  offspring_data = add_districts(offspring_data, cell_length = 1)
  offspring_data = offspring_data[offspring_data$district != 0,]
  
  # Add offspring events to data
  data = rbind(data, offspring_data)
  n_background[length(n_background) + 1] = nrow(data)
  
  # Set background data to offspring data for next iteration
  background_data = offspring_data
  background_data = background_data[background_data$t <= time_steps,]
}

# Remove out of time bounds data
data = data[data$t <= time_steps,]

# Print number of crimes
cat(paste0('Number of crimes: ', nrow(data), '\n'))
cat(paste0('Background crimes: ', nrow(data[data$parent_index == 0,]), '\n'))
cat(paste0('Aftershock crimes: ', nrow(data[data$parent_index != 0,]), '\n'))
cat(paste0('Maximal time step: ', max(data$t), '\n'))

# Plot superset crime in first time steps colored by district
# ggplot() + 
#   geom_bogota() + 
#   theme_bogota() + 
#   add_grid() +
#   geom_point(data = data[data$t <= 30,], aes(x = x, y = y, color = district))


t11 = Sys.time()
t_del = t11 - t00
cat(paste0('Superset aftershock done. Took ', t_del, ' ', units(t_del), '\n'))


###############
# Victimization data (thinned version of superset data)
###############

# Add thinning rates for each districts that give us the desired victimization in each district
n_semesters = time_steps/(365/2)
meta$keep_rate = 0
cat('Compute keep rates according to victimization...\n')
for (distr in as.character(meta$District)){
  desired_n_crimes = meta[meta$District == distr, 'n_crimes'] * n_semesters
  sampled_n_crimes = nrow(data[data$district == distr,])
  if (desired_n_crimes > sampled_n_crimes){
    stop(paste0('Not enough crimes for ', distr, '.\nSupersample: ', sampled_n_crimes, '\nDesired: ', desired_n_crimes))
  } else {
    cat('For ', distr, ': sampled_n = ', sampled_n_crimes, ', desired_n = ', desired_n_crimes, '\n')
  }
  meta[meta$District == distr, 'keep_rate'] = desired_n_crimes/sampled_n_crimes
}

# Thin accordingly
thin_data = data.frame(matrix(ncol = 5, nrow = 0))
for (distr in as.character(meta$District)){
  keep_rate = meta[meta$District == distr, 'keep_rate']
  rel_data = data[data$district == distr,]
  rownames(rel_data) = NULL
  bin = rbinom(n = nrow(rel_data), size = 1, prob = keep_rate)
  thin_data = rbind(thin_data, rel_data[which(bin == 1),])
}
rownames(thin_data) = NULL

# Plot thinned version of data in first 3 months
# ggplot() + 
#   geom_bogota() + 
#   geom_point(data = thin_data[thin_data$t < 60,], aes(x = x, y = y, color = district), alpha = 0.5) + 
#   ggtitle('Vicitmization data for 60 days')

# Save data
write.csv(thin_data, file = paste0('data/real_bogota_victim_data_v', version, '_ts', time_steps, '.csv'), row.names = FALSE)
write.csv(data, file = paste0('data/real_bogota_superset_data_v', version, '_ts', time_steps, '.csv'), row.names = FALSE)

t1 = Sys.time()
t_del = t1 - t0
cat(paste0('Done. Took ', t_del, ' ', units(t_del), '\n'))

