#############
# Simulation: Fit SEPP to Bogota data with and without thinning
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
library(plyr)
library(scales)
library(future.apply)

set.seed(0)

source('utils/model_utils.R')
source('utils/real_bogota_utils.R')

very_beginning = Sys.time()


############
# Adjust parameters
############

# Data parameters
time_steps = 6 * 365
version = 1 
data_path = paste0('data/real_bogota_victim_data_v', version, '_ts', time_steps, '.csv') # This is already thinned according to victimization rates

# Prefix for data folder 
output_folder = paste0('version_', version)

# Time steps for training
ignore_first_t = 500 
training_time_cutoff = 4 * 500

# Intial parameter guesses
init_mu_bar = 50
init_theta = 10
init_omega = 1
init_sigma_x = 0.5
init_sigma_y = 0.5

# Speed up -- truncation
crimes_per_day = 125
p_cut = 0.05
true_omega = 0.2
# Roughly: with prob at least 1-p_cut all offspring of crime in the training data occur within N timesteps
N = -log(1-((1-p_cut)^(1/(time_steps*crimes_per_day))))/true_omega
N = round(N)
cat(paste0('We assume offspring crimes occur within ', N, ' timesteps.\n'))
true_sigma = 0.1 # Assumes true_sigma_x == true_sigma_y
p_star = 0.5 * (1 - p_cut)^(1/(2*time_steps*crimes_per_day)) + 0.5 
eps = qnorm(p_star, mean = 0, sd = true_sigma)
cat(paste0('We assume all offspring occurs within ', eps, ' in x and y direction.\n'))
## The other way around for sanity
#pp = (pnorm(eps, mean = 0, sd = true_sigma) - pnorm(-eps, mean = 0, sd = true_sigma))^(2*crimes_per_day*time_steps)
#cat(paste0('Probability that all offspring lie within ', eps, ' of parent in both x an y is ', pp, '\n'))


############
# Meta data
############

# Read meta info: district, population, victimization (per semester), report rate (per semester)
meta_path = 'metadata/bogota_victimization.csv'
meta = read.csv(meta_path)

# Scale population down
meta$population_scaled = meta$Population / 40

# Number of crimes per semester (assume victimization equals exactly one crime)
meta$n_crimes = round(meta$population_scaled * meta$Victimization)

# Add thinning rate
meta$thinning_rate = 1 - meta$Percent_reported


############
# Read data
############

# Read data and restrict to training period 
data = read.csv(data_path)
data = data[,c('id', 'x', 'y', 't', 'parent_index', 'district')]
data = data[(data$t > ignore_first_t) & (data$t <= training_time_cutoff),]
rownames(data) = NULL
true_data = data
# hist(data$t)

# Create output folder
dir.create(paste0('output/', output_folder, '/'), showWarnings = FALSE)


############
# Start iteration 
############

for (thinning in c(TRUE, FALSE)) {
  
  cat(paste0('\nStart thinning = ', thinning, '\n'))
  
  # Reset data
  data = true_data
  
  # Path variables
  parameter_path = paste0('output/', output_folder, '/learned_parameters/learned_parameters_thinning_', thinning, '.csv')
  dir.create(paste0('output/', output_folder, '/learned_parameters/'), showWarnings = FALSE)
  
  training_data_path = paste0('output/', output_folder, '/training_data/training_data_thinning_', thinning, '.csv')
  training_data_plot_path = paste0('output/', output_folder, '/training_data/training_data_thinning_', thinning, '.png')
  dir.create(paste0('output/', output_folder, '/training_data/'), showWarnings = FALSE)
  
  
  #######################################################################
  # Add thinning
  #######################################################################
  
  if (thinning) {
    
    thin_data = data.frame(matrix(ncol = 6, nrow = 0))
    
    for (distr in as.character(meta$District)) {
      
      reportance_rate = meta[meta$District == distr, 'Percent_reported']
      rel_data = true_data[true_data$district == distr,]
      rownames(rel_data) = NULL
      bin = rbinom(n = nrow(rel_data), size = 1, prob = reportance_rate)
      thin_data = rbind(thin_data, rel_data[which(bin == 1),])
      cat(paste0('Number of unreported crimes in ', distr, ' : ', length(bin) - sum(bin), '/', length(bin), ' (', round(100*(length(bin) - sum(bin))/length(bin), digits = 2), '%)\n'))
      
    }
    
    # Set data to thinned data
    data = thin_data
    rownames(data) = NULL
  } 
  
  # Sanity check
  if (nrow(data) == 0){
    stop('Thinned data empty! This will throw error in evaluation script.')
  }
  
  # Save training data
  write.csv(data, file = training_data_path, row.names = FALSE)
  
  # Plot training
#  title = paste0('Training data thinning = ', thinning, ', t = [', ignore_first_t, ', ', training_time_cutoff,']')
#  training_data_plot = ggplot() + 
#    geom_bogota() +
#    theme_bogota() + 
#    geom_point(data = data, aes(x = x, y = y, color = district), alpha = 0.1) + 
#    ggtitle(title)
#  #training_data_plot
#  ggsave(filename = training_data_plot_path,
#         plot = training_data_plot,
#         width = 10,
#         height = 8,
#         dpi = 300,
#         units = "in")
  
  
  #######################################################################
  # Fit model with EM
  #######################################################################
  
  parameters = real_bogota_seppEM(data = data,
                      init_parameters = c(init_mu_bar, init_theta, init_omega, init_sigma_x, init_sigma_y),
                      cell_length = 1, 
                      max_iter_EM = 100,
                      tol_EM = 1e-5,
                      lower_EM = rep(0.01, 5),
                      upper_EM = c(100, 50, 10, 10, 10),
                      t_cutoff_EM = N,
                      xy_cutoff_EM = eps,
                      background_sd = 15)
  
  
  # Save learned parameters
  write.csv(parameters, file = parameter_path)
  
} #thinning = TRUE FALSE

very_end = Sys.time()
time_delta = very_end - very_beginning
cat(paste0('\n\nRun time total: ', time_delta, ' ', units(time_delta), '\n'))
