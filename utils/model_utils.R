########################################################
# Self-Exciting Spatio-Temporal Point Process Modeling
# Nil-Jana Akpinar, 2020
########################################################

library(pryr)
library(plyr)
library(assertthat)


############
# Helper functions
############

# Returns list of deterministic centers
get_centers = function(data_set){
  
  if (data_set == 'realistic_bogota') {
    
    list(c(6,20),
         c(-6,20),
         c(6,-20),
         c(-6,-20),
         c(6,-10),
         c(6,10),
         c(-6,-10),
         c(-6,10),
         c(-6,0),
         c(6,0),
         c(-6,-30),
         c(-6,30),
         c(6,30),
         c(6,-30))
    
  } else {
    stop('Invalid argument for get_centers!')
  }
}


############
# Define background intensity
############

# Specify unimodal background intensity
unimodal_background_intensity = function(x, y, center, mu_bar, sd = 4.5) {
  (mu_bar/(2*pi*sd**2))*exp(-(x-center[1])**2/(2*sd**2))*exp(-(y-center[2])**2/(2*sd**2))
}

# Specify unimodal background intensity, variable sigma
unimodal_background_intensity_variable_sigma = function(x, y, center, mu_bar, sigma) {
  (mu_bar/(2*pi*sigma**2))*exp(-(x-center[1])**2/(2*sigma**2))*exp(-(y-center[2])**2/(2*sigma**2))
}

# Returns unimodal background intensity function only dependent on (x,y)
get_unimodal_background_intensity_function = function(x, y, center, mu_bar){
  res_function = function(x,y) pryr::partial(unimodal_background_intensity, center = center, mu_bar = mu_bar)(x,y)
  res_function
}

# Returns unimodal background intensity function only dependent on (x,y), variable sigma
get_unimodal_background_intensity_function_variable_sigma = function(x, y, center, mu_bar, sigma){
  res_function = function(x,y) pryr::partial(unimodal_background_intensity_variable_sigma, center = center, mu_bar = mu_bar, sigma = sigma)(x,y)
  res_function
}

# Returns multimodal background intensity function only dependent on (x,y), need 5 centers as list 
get_multimodal_background_intensity_function = function(x, y, centers, total_mu_bar){
  
  # Centers can be found with get_centers()
  assert_that(is.list(centers))
  
  # Assume each of the parts equal
  mu_partial = total_mu_bar / length(centers)
  
  for (ii in 1:length(centers)){
    local({
      idx = ii
      assign(paste0('background_intensity_part', idx), get_unimodal_background_intensity_function(x, y, centers[[idx]], mu_partial), pos = .GlobalEnv)
    })
  }
  
  res_function = function(x,y){
      res = 0
      for (ii in 1:length(centers)){
        res = res + get(paste0('background_intensity_part', ii))(x,y)
      }
      res
    }

  # Return single background intensity function
  res_function
}


########################################################
# Define aftershock intensity
########################################################

# Aftershock intensity function, centered around (x_center, y_center, t_center)
aftershock_intensity = function(x, y, t, x_center, y_center, t_center, sigma_x, sigma_y, theta, omega){
  theta*omega*exp(-omega*(t-t_center))*exp(-(x-x_center)**2/(2*sigma_x**2))*exp(-(y-y_center)**2/(2*sigma_y**2))
}

# Get aftershock intensity function only dependent on (x,y,t)
get_aftershock_intensity = function(x, y, t, x_center, y_center, t_center, sigma_x, sigma_y, theta, omega){
  res_function = function(x,y,t) pryr::partial(aftershock_intensity,
                                               x_center = x_center,
                                               y_center = y_center,
                                               t_center = t_center,
                                               sigma_x = sigma_x,
                                               sigma_y = sigma_y,
                                               theta = theta,
                                               omega = omega)(x,y,t)
  res_function
}


########################################################
# Define total intensity function -- lambda
########################################################

# Returns total intensity function only dependent on (x, y, t) OR dependent on (x, y, t, data)
get_intensity_function = function(x, y, t, centers, total_mu_bar, sigma_x, sigma_y, theta, omega, data = NULL){
  
  # Background intensity function, multimodal
  background_int = get_multimodal_background_intensity_function(x, y, centers, total_mu_bar)
  
  # Aftershock intensity, uncentered
  aftershock_int = get_aftershock_intensity(x, y, t, x_center = 0, y_center = 0, t_center = 0, sigma_x, sigma_y, theta, omega)
  
  if (empty(data)){
    # if: function should take data as argument
    res_function = function(x, y, t, data){
      rel_data = data[data$t < t,]
      res1 = sum(aftershock_int( x - rel_data$x, y - rel_data$y, t - rel_data$t))
      res2 = background_int(x, y)
      res1 + res2
    }
  } else {
    # else: function should already incorporate the data
    # Sum over previous events
    res_function = function(x, y, t){
      rel_data = data[data$t < t,]
      res1 = sum(aftershock_int( x - rel_data$x, y - rel_data$y, t - rel_data$t))
      res2 = background_int(x, y)
      res1 + res2
    }
  }
  
  # Return function
  res_function
}