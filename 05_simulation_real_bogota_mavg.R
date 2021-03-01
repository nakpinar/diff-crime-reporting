#############
# Simulation: Fit MAVG to Bogota data with and without thinning
# Nil-Jana Akpinar, 2020
#############

library(class)
library(pryr)
library(ggplot2)
library(FNN)
library(cubature)
library(VGAM)
library(reshape2)
library(svMisc)
library(Matrix)
library(plyr)
library(future.apply)

set.seed(0)

source('utils/real_bogota_utils.R')

very_beginning = Sys.time()

#############
# Adjust parameters
#############

version = 1
output_folder = paste0('version_', version)
dir.create(paste0('output/', output_folder, '/mavg_parameters/'), showWarnings = FALSE)


#############
# Start iteration
#############

for (thinning in c(TRUE, FALSE)) {
  
  cat(paste0('\nStart thinning = ', thinning, '\n'))
  
  # Read training data (from SEPP fitting)
  training_data_path = paste0('output/', output_folder, '/training_data/training_data_thinning_', thinning, '.csv')
  data = read.csv(training_data_path)
  
  ## Get alphas to search over
  # Bandwidth to omega
  Tt = max(data$t) - min(data$t) + 1
  x_T = max(data$t)
  x_1 = min(data$t)
  h_to_w = function(h) exp(-h * (Tt - 1) / (x_T - x_1))
  
  days_to_search = 1:50
  h = 1/days_to_search
  ws = h_to_w(h)
  alphas = 1 - ws
  #plot(alphas, days_to_search)
  
  # Fit model
  res = real_bogota_exp_mavg(data = data,
                       alphas = alphas,
                       cell_length = 1)
  
  losses = data.frame(alpha = alphas, loss = res[[2]])
  parameter = data.frame(alpha = res[[1]])
  
  # Save results
  parameter_path = paste0('output/', output_folder, '/mavg_parameters/learned_parameter_thinning_', thinning, '.csv')
  losses_path = paste0('output/', output_folder, '/mavg_parameters/losses_thinning_', thinning, '.csv')
  
  write.csv(parameter, file = parameter_path, row.names = FALSE)
  write.csv(losses, file = losses_path, row.names = FALSE)
  
  #alpha_to_days = function(alpha) (1 / (-(x_T - x_1) / (Tt - 1) * log(1 - alpha))):
  #r = alpha_to_days(res[[1]])  
  #cat('Best paramater =', parameter, '\nCorresponds to', alpha_to_days(r), 'days')
  } 

very_end = Sys.time()
time_delta = very_end - very_beginning
cat(paste0('\n\nRun time total: ', time_delta, ' ', units(time_delta), '\n'))
