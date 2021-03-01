#######################################################################
# Auxiliary functions for Bogota simulation; realistic map
# Nil-Jana Akpinar, 2020
#######################################################################

library(ggplot2)
library(assertthat)
library(rgeos)
library(grid)
library(VGAM)
library(fastDummies)
library(forecast)
library(dplyr)
library(sf)


############
# Meta data 
############

# Read Bogota meta files
meta_path = 'metadata/bogota_victimization.csv'
meta = read.csv(meta_path)
bogota_shp = 'metadata/bogota.shp'
bogota = st_read(bogota_shp, stringsAsFactors = FALSE)

# Read bogota_mask dependent on cell_length
get_bogota_mask = function(cell_length){
  
  bogota_mask = read.csv(paste0('metadata/bogota_mask_', cell_length, '.csv'))
  bogota_mask$cell_id = 1:nrow(bogota_mask)
  bogota_mask
  
}

# Add which hot spot cell (i.e. mask_id) each crimes lies in
add_mask_cells = function(data, cell_length){
  
  mask = get_bogota_mask(cell_length = cell_length)
  
  data$x_ceil = ceiling(data$x)
  data$y_ceil = ceiling(data$y)
  mask$x_ceil = mask$x
  mask$y_ceil = mask$y
  res = merge(x = data, y = mask[, c('x_ceil', 'y_ceil', 'district', 'cell_id')], all.x = TRUE)
  res[, - which(colnames(res) %in% c('x_ceil', 'y_ceil'))]
  
}

# Return xmin,xmax,ymin,ymax for broad Bogota plotting area
get_bogota_window = function(data = bogota){
  
  # Get maxmin boundaries of districts
  bbox_list = lapply(st_geometry(data), st_bbox)
  maxmin = as.data.frame(matrix(unlist(bbox_list), nrow = nrow(data), byrow = TRUE))
  colnames(maxmin) = c('xmin', 'ymin', 'xmax', 'ymax')
  maxmin$locality = bogota$LocNombre
  
  # Full maximum range of the city
  xmin = min(maxmin$xmin)
  xmax = max(maxmin$xmax)
  ymin = min(maxmin$ymin)
  ymax = max(maxmin$ymax)
  
  list(c(xmin, xmax), c(ymin, ymax))
}

# Takes data frame and adds districts column
add_districts = function(data, cell_length){
  
  mask = get_bogota_mask(cell_length = cell_length)
  
  data$x_ceil = ceiling(data$x)
  data$y_ceil = ceiling(data$y)
  mask$x_ceil = mask$x
  mask$y_ceil = mask$y
  res = merge(x = data, y = mask[, c('x_ceil', 'y_ceil', 'district')], all.x = TRUE)
  res[, - which(colnames(res) %in% c('x_ceil', 'y_ceil'))]

}

# Adds true background crime, i.e. absolute cluster to data
add_true_backgrounds = function(data){
  
f =  function(i){
  row = data[i,]
  print(i) #History
  print(row)
  if (row$parent_index == 0) { 
    return(row$id)
  } else {
    f(which(data$id == row$parent_index))
  }
}

plan(multiprocess)
YY = future_lapply(1:nrow(data), f)

data$true_background = unlist(YY)
data
}


############
# ggplot2 components for Bogota plots
############

# Gets n ggplot colors
get_colors = function(n){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

theme_bogota = function(){
  theme(panel.grid.major = element_line(colour = "transparent"), panel.background = element_rect(fill = NA, color = 'black'))
}

geom_bogota = function(data = bogota, fill_district = FALSE){
  if (fill_district){
    geom_sf(data = data, aes(fill = LocNombre))
  } else {
    geom_sf(data = data, fill = alpha('white', 0)) 
  }
}

add_grid = function(cell_length = 1, fill_district = FALSE, alpha = 0.6, add_gridlines = FALSE){
  
  mask = get_bogota_mask(cell_length = 1)
  
  if (fill_district){
    if (add_gridlines){
      list(geom_raster(data = mask, aes(x = x, y = y, fill = district), alpha = alpha, hjust = 0, vjust = 0),
           geom_hline(yintercept = mask$y, size = 0.3),
           geom_vline(xintercept = mask$x, size = 0.3),
           scale_y_continuous(limits = c(min(mask$y) - 1, max(mask$y) + 1), expand = c(0, 0)),
           scale_x_continuous(limits = c(min(mask$x) - 1, max(mask$x) + 1), expand = c(0, 0)))
    } else {
      list(geom_raster(data = mask, aes(x = x, y = y, fill = district), alpha = alpha, hjust = 0, vjust = 0),
           scale_y_continuous(limits = c(min(mask$y) - 1, max(mask$y) + 1), expand = c(0, 0)),
           scale_x_continuous(limits = c(min(mask$x) - 1, max(mask$x) + 1), expand = c(0, 0)))
    }
  } else {
    list(geom_hline(yintercept = mask$y, size = 0.3),
         geom_vline(xintercept = mask$x, size = 0.3),
         scale_y_continuous(limits = c(min(mask$y) - 1, max(mask$y) + 1), expand = c(0, 0)),
         scale_x_continuous(limits = c(min(mask$x) - 1, max(mask$x) + 1), expand = c(0, 0)))
  }
  
}


############
# Expectation Maximization for SEPP in Bogota bounds
############

# Returns polygon of the outer grid lines (i.e. Bogota outline)
get_outer_grid_polygon = function(cell_length){
  
  bogota_mask = get_bogota_mask(cell_length = cell_length)
  
  grid_polygon_points = data.frame(x = c(-13, -11, -11, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 3, 3, 4, 4, 5, 5, 7, 7, 8, 8, 9, 9, 11, 11, 12, 12, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9,
                   11, 11, 12, 12, 11, 11, 10, 10, 11, 11, 10, 10, 11, 11, 12, 12, 11, 11, 7, 7, 2, 2, 0, 0, -1, -1, -2, -2, -3, -3, -4, -4, -6, -6, -8, -8, -7, -7, -8, -8, -7, -7,
                   -8, -8, -11, -11, -12, -12, -13, -13, -12, -12, -11, -11, -9, -9, -8, -8, -9, -9, -8, -8, -9, -9, -10, -10, -11, -11, -12, -12, -13, -13, -12, -12, -13, -13),
             y = c(-31,-31, -30, -30, -29, -29, -28, -28, -27, -27, -26, -26, -24, -24, -19, -19, -13, -13, -11, -11, -8, -8, -7, -7, -5, -5, -3, -3, -1, -1, 1, 1, 4, 4, 7, 7, 9, 9, 8, 8, 9, 9, 10, 10, 12,
                   12, 13, 13, 14, 14, 15, 15, 17, 17, 23, 23, 26, 26, 27, 27, 28, 28, 30, 30, 31, 31, 27, 27, 26, 26, 23, 23, 21, 21, 20, 20, 19, 19, 17, 17, 15, 15, 14, 14, 13, 13, 12,
                   12, 11, 11, 10, 10, 9, 9, 7, 7, 6, 6, 5, 5, 4, 4, -2, -2, -4, -4, -8, -8, -12, -12, -14, -14, -16, -16, -21, -21, -26, -26, -30, -30, -31))
  polygon_string = paste0('LINESTRING(')
  for (row in 1:nrow(grid_polygon_points)){
    if (row != nrow(grid_polygon_points)){
      polygon_string = paste0(polygon_string, grid_polygon_points[row,'x'], ' ', grid_polygon_points[row,'y'], ',')
    } else {
      polygon_string = paste0(polygon_string, grid_polygon_points[row,'x'], ' ', grid_polygon_points[row,'y'], ')')
    }
  }
  
  polygon = readWKT(polygon_string)
  # plot(polygon)
  polygon
}

# Adds the minimum Euclidean distance to Bogota edge (as defined by grid) to each point in data
add_dist_edge = function(data, cell_length){
  
  t0 = Sys.time()
  cat('Compute distances to grid edges...\n')
  
  grid_edge = get_outer_grid_polygon(cell_length = cell_length)
  spatial_points = SpatialPoints(data[, c('x', 'y')])
  
  f = function(row){
    gDistance(spgeom1 = spatial_points[row], spgeom2 = grid_edge)
  }
  
  res = unlist(lapply(1:length(spatial_points), f))
  data$dist_grid_edge = res
  
  t1 = Sys.time()
  time_delta = t1 - t0
  cat(paste0('Done. Took: ', time_delta, ' ', units(time_delta), '\n'))
  
  data
}

# Get int_{bogota bounds} (x in N(x_center, sigma_x))(y in N(y_center, sigma_y)) d(x,y)
# Mask points give the upper right corner of cell; integral is over mask cells not over polygon.
# VERSION 1: for one center point center = c(x,y)
get_bogota_normal_int = function(center, sds, cell_length, mask){
  
  bogota_mask_help = mask[mask$district != 0,]

  xseq = unique(bogota_mask_help$x)
  xseq_ints = pnorm(xseq, mean = center[1], sd = sds[1]) - pnorm(xseq - cell_length, mean = center[1], sd = sds[1])
  x_ints = as.data.frame(cbind(xseq, xseq_ints))
  colnames(x_ints) = c('x', 'x_ints')
  bogota_mask_help = merge(bogota_mask_help, x_ints, all.x = TRUE)

  yseq = unique(bogota_mask_help$y)
  yseq_ints = pnorm(yseq, mean = center[2], sd = sds[2]) - pnorm(yseq - cell_length, mean = center[2], sd = sds[2])
  y_ints = as.data.frame(cbind(yseq, yseq_ints))
  colnames(y_ints) = c('y', 'y_ints')
  bogota_mask_help = merge(bogota_mask_help, y_ints, all.x = TRUE)

  bogota_mask_help$prod_int = bogota_mask_help$x_ints * bogota_mask_help$y_ints
  sum(bogota_mask_help$prod_int)

}

# VERSION 2: for a df of center points centers = df with cols x,y
get_bogota_normal_int_data = function(centers, sds, cell_length){ # centers is df with x y cols
  
  mask = get_bogota_mask(cell_length = cell_length)
  bogota_mask_help = mask[mask$district != 0,]

  xseq = unique(bogota_mask_help$x)
  xseq_mat = matrix(rep(xseq, each = nrow(centers)), nrow = nrow(centers))
  res_x = 0.5 * (1 + erf((xseq_mat - centers$x)/(sds[1] * sqrt(2)))) - 0.5 * (1 + erf((xseq_mat - cell_length - centers$x)/(sds[1] * sqrt(2))))
  
  yseq = unique(bogota_mask_help$y)
  yseq_mat = matrix(rep(yseq, each = nrow(centers)), nrow = nrow(centers))
  res_y = 0.5 * (1 + erf((yseq_mat - centers$y)/(sds[2] * sqrt(2)))) - 0.5 * (1 + erf((yseq_mat - cell_length - centers$y)/(sds[2] * sqrt(2))))
  
  res = matrix(0, nrow = nrow(centers), ncol = nrow(bogota_mask_help))
  for (i in 1:nrow(bogota_mask_help)){
    x = bogota_mask_help[i,'x']
    y = bogota_mask_help[i,'y']
    res[,i] = res_x[,which(xseq == x)] * res_y[,which(yseq == y)]
  }
  
  rowSums(res)
}

# VERSION 3: for a df of center points data = df with cols x,y. Check if points are at the edge and worth considering first.
get_bogota_normal_int_data_edge_check = function(data, sds, cell_length){

  assert_that('x' %in% colnames(data))
  assert_that('y' %in% colnames(data))
  
  if (!('dist_grid_edge' %in% colnames(data))){
    data = add_dist_edge(data = data, cell_length = cell_length)
  }
  
  # We consider the integrals over the largest square in the dist to edge circle.
  # If that is too small, the integral over all Bogota needs to be computed, otherwise it's 1.
  data$half_sl = data$dist_grid_edge / sqrt(2) # Half the side length of the square
  x_square_ints = pnorm(data$half_sl, mean = 0, sd = sds[1]) - pnorm(-data$half_sl, mean = 0, sd = sds[1])
  y_square_ints = pnorm(data$half_sl, mean = 0, sd = sds[2]) - pnorm(-data$half_sl, mean = 0, sd = sds[2])
  square_ints = x_square_ints * y_square_ints
  # Index of the points that do not have integral 1 in the small square; if too slow we can adjust threshold here (with 1 roughly 15% of centers left).
  idx = which(square_ints < 0.99)
  
  # Default for aftershock normal int is 1
  data$aftershock_int = 1
  
  if (length(idx) > 0){
    bogota_mask = get_bogota_mask(cell_length = cell_length)
    
    data[idx, 'aftershock_int'] = get_bogota_normal_int_data(centers = data[idx,c('x','y')],
                               sds = sds,
                               cell_length = 1)
  }
    
  data$aftershock_int
}

# Fits SEPP with unimodal background with EM
real_bogota_seppEM = function(data, # Data frame with (at least) columns t,x,y
                  init_parameters, # Initial paramaters; c(init_mu_bar, init_theta, init_omega, init_sigma_x, init_sigma_y)
                  cell_length = 1,
                  max_iter_EM = 100,
                  tol_EM = 0.001,
                  lower_EM = rep(0.01, 5), # Box constraints for optimization of the parameters
                  upper_EM = c(100, 50, 10, 10, 10),
                  t_cutoff_EM = FALSE, # Time and space cut-offs for estimation of the branching structure; values have to be computed externally and given here 
                  xy_cutoff_EM = FALSE,
                  background_sd = 15 # Std of background model (for realistic Bogota map 15 is appropriate)
){
  
  # Input assertions
  assert_that(is.data.frame(data))
  assert_that('x' %in% colnames(data))
  assert_that('y' %in% colnames(data))
  assert_that('t' %in% colnames(data))
  assert_that(is.vector(init_parameters))
  assert_that(length(init_parameters) == 5)
  assert_that(all(lower_EM < upper_EM))
  assert_that(all(lower_EM < init_parameters))
  assert_that(all(init_parameters < upper_EM))
  assert_that(is.number(t_cutoff_EM) | (t_cutoff_EM == FALSE))
  assert_that(is.number(xy_cutoff_EM) | (xy_cutoff_EM == FALSE))
  
  # Get mask 
  bogota_mask = get_bogota_mask(cell_length = cell_length)
  
  # Add distance to edge 
  data = add_dist_edge(data, cell_length = cell_length)
  
  # Shift t (such that we don't have to worry about the ignore_first_t)
  if (min(data$t) != 0){
    data$t = data$t - min(data$t)
  }
  
  # Get model with unimodal background around (0,0) and uncentered aftershock intensity
  center = c(0,0)
  model_background_intensity = function(x, y, mu_bar) pryr::partial(unimodal_background_intensity, center = center, sd = background_sd)(x, y, mu_bar)
  model_aftershock_intensity = function(x, y, t, sigma_x, sigma_y, theta, omega) pryr::partial(aftershock_intensity, x_center = 0, y_center = 0, t_center = 0)(x, y, t, sigma_x, sigma_y, theta, omega)
  
  # Set initial values for parameter estimates and tracked variabels
  mu_bars = c(init_parameters[1]) 
  thetas = c(init_parameters[2])
  omegas = c(init_parameters[3]) 
  sigma_xs = c(init_parameters[4]) 
  sigma_ys = c(init_parameters[5]) 
  
  Q_ttm1= -1
  Q_tm1tm1 = -1
  vec_sum_pui = -1
  vec_sum_pui_t_diff = -1
  vec_sum_pui_x_diff_sq = -1
  vec_sum_pui_y_diff_sq = -1
  vec_sum_parwise_g = -1
  vec_sum_pui0 = -1
  
  iter = 1
  data = data[order(data$t),]
  
  # Pre-compute int_{bogota bounds} (x in N(0,15))(y in N(0,15)) d(x,y)
  bogota_bckg_int = get_bogota_normal_int(center = c(0,0), sds = c(background_sd,background_sd), cell_length = 1, mask = bogota_mask)
  
  while(iter < max_iter_EM){
    
    t0iter = Sys.time()
    cat(paste0('\nEM iteration: ', iter, '/max', max_iter_EM, '\n'))
    
    # Set parameters for this iteration
    mu_bar = mu_bars[iter]
    theta = thetas[iter]
    omega = omegas[iter]
    sigma_x = sigma_xs[iter]
    sigma_y = sigma_ys[iter]
    
    cat('Parameters for this iteration: mu_bar = ', mu_bar, ', theta = ', theta, ', omega = ', omega, ', sigma_x = ', sigma_x, ', sigma_y = ', sigma_y, '\n')
    
    # Stop if parameter estimates converged
    if(iter > 1){
      if((abs(mu_bar - mu_bars[iter - 1]) < tol_EM) &
         (abs(theta - thetas[iter - 1]) < tol_EM) &
         (abs(omega - omegas[iter - 1]) < tol_EM) &
         (abs(sigma_x - sigma_xs[iter - 1]) < tol_EM) &
         (abs(sigma_y - sigma_ys[iter - 1]) < tol_EM)){
        cat(paste0('\nIterations: ',iter))
        break
      }
    }
    
    # Compute latent variables based on above parameters
    ## Latent pui and pui-based variables
    cat('Compute pairwise triggering probabilities P(u_i=j)...\n')
    t0 = Sys.time()
    
    plan(multiprocess)
    
    YY = future_lapply(1:nrow(data), function(ii){
      
      ii_t = data[ii,'t']
      ii_id = data[ii,'id']
      ii_x = data[ii, 'x']
      ii_y = data[ii, 'y']
      
      if ((t_cutoff_EM != FALSE) & (xy_cutoff_EM != FALSE)){
        rel_data = data[(data$t < ii_t) & (ii_t - data$t < t_cutoff_EM) & (abs(ii_x - data$x) < xy_cutoff_EM) & (abs(ii_y - data$y) < xy_cutoff_EM),]
      } else if (t_cutoff_EM != FALSE) {
        rel_data = data[(data$t < ii_t) & (ii_t - data$t < t_cutoff_EM),]
      } else if (xy_cutoff_EM != FALSE){
        rel_data = data[(data$t < ii_t) & (abs(ii_x - data$x) < xy_cutoff_EM) & (abs(ii_y - data$y) < xy_cutoff_EM),]
      } else {
        rel_data = data[(data$t < ii_t),]
      }
      
      if (nrow(rel_data) > 0){
        ints = model_aftershock_intensity(x = data$x[ii] - rel_data$x,
                                          y = data$y[ii] - rel_data$y,
                                          t = data$t[ii] - rel_data$t,
                                          sigma_x = sigma_x,
                                          sigma_y = sigma_y,
                                          theta = theta,
                                          omega = omega)
        sum_ints = sum(ints)
        
        # pui 
        mu = model_background_intensity(ii_x, ii_y, mu_bar)
        lam = mu + sum_ints
        pui = ints / lam 
        sum_pui = sum(pui)
        
        # pui * t_diff
        t_diff = ii_t - rel_data$t
        sum_pui_t_diff = sum(t_diff * pui)
        
        # pui * x_diff**2
        x_diff = ii_x - rel_data$x
        sum_pui_x_diff_sq = sum(x_diff ** 2 * pui)
        
        # pui * y_diff**2
        y_diff = ii_y - rel_data$y
        sum_pui_y_diff_sq = sum(y_diff ** 2 * pui)
      } else {
        sum_ints = 0 
        sum_pui = 0
        sum_pui_t_diff = 0
        sum_pui_x_diff_sq = 0
        sum_pui_y_diff_sq = 0
      }
      
      c(sum_ints, sum_pui, sum_pui_t_diff, sum_pui_x_diff_sq, sum_pui_y_diff_sq)
    })
    
    YY_sum = Reduce('+', YY) 
    sum_pui = YY_sum[2]
    sum_pui_t_diff = YY_sum[3]
    sum_pui_x_diff_sq = YY_sum[4]
    sum_pui_y_diff_sq = YY_sum[5]
    sum_pairwise_g = sapply(YY, "[[", 1)
    
    ## Latent pui0
    mu = model_background_intensity(data$x, data$y, mu_bar)
    lambda_mat = mu + sum_pairwise_g
    pui0 = mu/lambda_mat
    sum_pui0 = sum(pui0)
    
    # Add to tracking 
    vec_sum_pui[iter + 1] = sum_pui
    vec_sum_pui_t_diff[iter + 1] = sum_pui_t_diff
    vec_sum_pui_x_diff_sq[iter + 1] = sum_pui_x_diff_sq
    vec_sum_pui_y_diff_sq[iter + 1] = sum_pui_y_diff_sq
    vec_sum_parwise_g[iter + 1] = sum(sum_pairwise_g)
    vec_sum_pui0[iter + 1] = sum_pui0
    
    t1 = Sys.time()
    t_del = t1 - t0
    cat(paste0('Done. Took ', t_del, ' ', units(t_del), '\n'))
    
    # Q function -- Eloglik complete with mu and the constant term
    Q = function(arg){
      mu_bar = arg[1]
      theta = arg[2]
      omega = arg[3]
      sigma_x = arg[4]
      sigma_y = arg[5]
      
      aftershock_int_colname = paste0('aftershock_integrals_', sigma_x, '_', sigma_y)
      data[,aftershock_int_colname] = get_bogota_normal_int_data_edge_check(data = data, sds = c(sigma_x, sigma_y), cell_length = cell_length)
      ints = data[, aftershock_int_colname]
  
      # Loglik
      one_two = log(theta) * sum_pui + log(omega) * sum_pui - omega * sum_pui_t_diff - (1/(2*sigma_x**2)) * sum_pui_x_diff_sq - (1/(2*sigma_y**2)) * sum_pui_y_diff_sq
      three = theta*2*pi*sigma_x*sigma_y*sum((1-exp(-omega*(max(data$t) - data$t))) *
                                                 ints)
      mu_term = log(mu_bar) * sum_pui0 - mu_bar * max(data$t) * bogota_bckg_int
      const =  sum_pui0 * log(2 * pi * background_sd ** 2) + sum(pui0 * (data$x^2 + data$y^2) / (2 * background_sd ** 2))
      res = one_two - three + mu_term - const
      res
    }
    
    # Optimization step
    t0 = Sys.time()
    Q_tm1tm1[iter + 1] = Q(c(mu_bar, theta, omega, sigma_x, sigma_y))
    cat('Optimizing...\n')
    opt = optim(par = c(mu_bar, theta, omega, sigma_x, sigma_y),
                fn = Q,
                control = list(fnscale = -1),
                method ='L-BFGS-B',
                lower = lower_EM,
                upper = upper_EM)
    par = opt$par
    t1 = Sys.time()
    t_del = t1 - t0
    cat(paste0('Done. Took ', t_del, ' ', units(t_del), '\n'))
    
    # Save new parameter estimates and ramaining tracking values
    mu_bars[iter + 1] = par[1]
    thetas[iter + 1] = par[2]
    omegas[iter + 1] = par[3]
    sigma_xs[iter + 1] = par[4]
    sigma_ys[iter + 1] = par[5]
    Q_ttm1[iter + 1] = Q(par)
    
    t1iter = Sys.time()
    t_del = t1iter - t0iter
    cat(paste0('Iteration step took ', t_del, ' ', units(t_del), '\n'))
    
    iter = iter + 1
  } # END while iter
  
  # Print learned parameters
  parameters = as.data.frame(cbind(mu_bars,
                                   thetas,
                                   omegas,
                                   sigma_xs,
                                   sigma_ys,
                                   Q_tm1tm1,
                                   Q_ttm1,
                                   vec_sum_pui,
                                   vec_sum_pui_t_diff,
                                   vec_sum_pui_x_diff_sq,
                                   vec_sum_pui_y_diff_sq,
                                   vec_sum_parwise_g,
                                   vec_sum_pui0))
  colnames(parameters) = c('mu_bar',
                           'theta',
                           'omega',
                           'sigma_x',
                           'sigma_y',
                           'Q(Theta_{t-1}, Theta_{t-1})',
                           'Q(Theta_{t}, Theta_{t-1})',
                           'sum_pui',
                           'sum_pui_t_diff',
                           'sum_pui_x_diff',
                           'sum_pui_y_diff',
                           'sum_pairwise_g',
                           'sum_pui0')
  
  parameters
}


############
# Within-cell moving average in Bogota bounds 
############

# Gaussian kernel function
gaussian_kernel = function(t) (1/sqrt(2*pi)) * exp(-t^2/2)

# Exponential kernel function
exp_kernel = function(t) ifelse(t <= 0, exp(t), 0)

# Fits exponential smoothing predictor for within cell MAVG
# alphas: search list of parameter values in [0,1], alpha = 1 - omega
real_bogota_exp_mavg = function(data, alphas, cell_length, district = FALSE){
  
  t_start = Sys.time()
  
  # Add cell_ids for each crime
  if (!('cell_id' %in% colnames(data))){
    data = add_mask_cells(data = data,
                          cell_length = cell_length)
  }
  
  # Fit MAVG model
  cat('Fit kernel MAVG model...\n')
  
  losses = c()
  
  bogota_mask = get_bogota_mask(cell_length = cell_length)
  bogota_mask_help = bogota_mask[bogota_mask$district != 0,]
  
  if (district != FALSE){
    assert_that(unique(data$district)[1] == district)
    bogota_mask_help = bogota_mask_help[bogota_mask_help$district == district,]
  }
  
  
  # Compute loss for each time_window
  for (alpha in alphas){
    
    t0 = Sys.time()
    cat('Compute loss for alpha = ', alpha, '...')
    
    # Number of crimes per day per cell (and add zero counts)
    data$count = 1
    crime_counts = aggregate(list(count = data$count), by = list(t = data$t, cell_id = data$cell_id), FUN = sum)
    all_combis = expand.grid(t = min(data$t):max(data$t), cell_id = bogota_mask_help$cell_id)
    crime_counts = merge(x = crime_counts, y = all_combis, by = c('cell_id', 't'), all.y = TRUE)
    crime_counts$count[is.na(crime_counts$count)] = 0
    
    # Sort by cell_id and then by t
    crime_counts = crime_counts[with(crime_counts, order(cell_id, t)), ]
    
    # Add cell_wise predictions with given alpha
    add_cell_wise_pred = function(x){
      x %>% mutate(pred = as.vector(ses(count, alpha = alpha, h = 1, initial = 'simple')$fitted))
    }
    
    crime_counts = as.data.frame(crime_counts %>%
                    group_split(cell_id) %>%
                    purrr::map_dfr(add_cell_wise_pred))
    
    # Compute MSE loss
    loss = sum((crime_counts$count - crime_counts$pred)^2)/nrow(crime_counts)
    cat(' Loss = ', loss, '\n')
    losses[length(losses) + 1] = loss
    
  }
  
  t_end = Sys.time()
  time_delta = t_end - t_start
  cat('MAVG fitting done. Took', time_delta, ' ', units(time_delta), '.\n')
  
  # Returns list with [[1]] = optimal time window, [[2]] = mean loss for each time window, [[3]] = loss over time steps for each time window
  list(alphas[which.min(losses)], losses)
}


############
# Integration for true intensity function
############

# Helper function for 01_true_integrals: Gets background integrals on grid with current value of mu.
# This is for multimodal background -- i.e. true integrals.
get_background_integrals = function(centers, cell_length, total_mu_bar, background_sd = 4.5){
  
  mask = get_bogota_mask(cell_length = cell_length)
  bogota_mask_help = mask[mask$district != 0,]
  
  xseq = unique(bogota_mask_help$x)
  xseq_mat = matrix(rep(xseq, each = nrow(centers)), nrow = nrow(centers))
  res_x = 0.5 * (1 + erf((xseq_mat - centers$x)/(background_sd * sqrt(2)))) - 0.5 * (1 + erf((xseq_mat - cell_length - centers$x)/(background_sd * sqrt(2))))
  
  yseq = unique(bogota_mask_help$y)
  yseq_mat = matrix(rep(yseq, each = nrow(centers)), nrow = nrow(centers))
  res_y = 0.5 * (1 + erf((yseq_mat - centers$y)/(background_sd * sqrt(2)))) - 0.5 * (1 + erf((yseq_mat - cell_length - centers$y)/(background_sd * sqrt(2))))
  
  res = matrix(0, nrow = nrow(centers), ncol = nrow(bogota_mask_help))
  for (i in 1:nrow(bogota_mask_help)){
    x = bogota_mask_help[i,'x']
    y = bogota_mask_help[i,'y']
    res[,i] = res_x[,which(xseq == x)] * res_y[,which(yseq == y)]
  }
  
  bogota_mask_help$superset_background_int = colSums(res) * total_mu_bar / nrow(centers)
  bogota_mask_help
}

# Helper function for 01_true_integrals: Gets sum(aftershock_int(relevant_crimes)) given a time step ts.
get_aftershock_integrals = function(centers, ts, sds, theta, omega, cell_length){
  
  mask = get_bogota_mask(cell_length = cell_length)
  bogota_mask_help = mask[mask$district != 0,]
  
  xseq = unique(bogota_mask_help$x)
  xseq_mat = matrix(rep(xseq, each = nrow(centers)), nrow = nrow(centers))
  res_x = 0.5 * (1 + erf((xseq_mat - centers$x)/(sds[1] * sqrt(2)))) - 0.5 * (1 + erf((xseq_mat - cell_length - centers$x)/(sds[1] * sqrt(2))))
  
  yseq = unique(bogota_mask_help$y)
  yseq_mat = matrix(rep(yseq, each = nrow(centers)), nrow = nrow(centers))
  res_y = 0.5 * (1 + erf((yseq_mat - centers$y)/(sds[2] * sqrt(2)))) - 0.5 * (1 + erf((yseq_mat - cell_length - centers$y)/(sds[2] * sqrt(2))))
  
  res = matrix(0, nrow = nrow(centers), ncol = nrow(bogota_mask_help))
  for (i in 1:nrow(bogota_mask_help)){
    x = bogota_mask_help[i,'x']
    y = bogota_mask_help[i,'y']
    res[,i] = res_x[,which(xseq == x)] * res_y[,which(yseq == y)] * theta * omega * exp(-omega * (ts - centers$t))
  }
  
  bogota_mask_help$superset_aftershock_int = colSums(res) * 2 * pi * sds[1] * sds[2]
  bogota_mask_help
}


############
# Ranking metric and hot spot plots
############

# Cumulative gain -- Sum of true integrated intensities in the top k predicted hotspots of modeled integrated intensities
cg_k = function(true_ints, estim_ints, k){
  
  assert_that(nrow(true_ints) == nrow(estim_ints))
  assert_that('victim_int' %in% colnames(true_ints))
  assert_that('int' %in% colnames(estim_ints))
  assert_that(k <= ceiling(0.2 * nrow(true_ints)))
  
  estim_rank = merge(x = estim_ints, y = true_ints[, c('x', 'y', 'district', 'victim_int')], all.x = TRUE)
  estim_rank = estim_rank[order(-estim_rank$int),]
  rownames(estim_rank) = NULL
  
  sum(estim_rank[1:k, 'victim_int'])
  
}

# Hot spot plots
plot_hotspots = function(true_ints, estim_ints, title, alpha = 0.4, k = 100, add_cg = TRUE, add_district_labels = FALSE){
  
  assert_that(nrow(true_ints) == nrow(estim_ints))
  assert_that('victim_int' %in% colnames(true_ints))
  assert_that('int' %in% colnames(estim_ints))
  assert_that(k <= ceiling(0.2 * nrow(true_ints)))
  
  # Add masks indicating top k values ( * alpha)
  estim_ints = estim_ints[order(-estim_ints$int),]
  rownames(estim_ints) = NULL
  estim_ints$topk = 0
  estim_ints$topk[1:k] = alpha
  
  true_ints = true_ints[order(-true_ints$victim_int),]
  rownames(true_ints) = NULL
  true_ints$topk = 0
  true_ints$topk[1:k] = alpha
  
  colors = get_colors(2)
  
  p  = ggplot() +
    geom_bogota() + 
    theme_bogota() +
    geom_raster(data = estim_ints, aes(x = x, y = y, fill = colors[1]), alpha = estim_ints$topk, hjust = 0, vjust = 0) + 
    geom_raster(data = true_ints, aes(x = x, y = y, fill = colors[2]), alpha = true_ints$topk, hjust = 0, vjust = 0) +
    scale_fill_discrete(name = '', labels = c('True', 'Predicted')) + 
    theme(legend.position=c(0.25, 0.9),
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10),
          legend.key.size = unit(1, "lines")) + 
    ggtitle(title) + 
    xlab('') + 
    ylab('')
  
  if (add_cg){
    
    cg = round(cg_k(true_ints = true_ints, estim_ints = estim_ints, k = k), 2)
    true_ints = true_ints[order(-true_ints$victim_int),]
    cg_star = round(sum(true_ints[1:k, 'victim_int']), 2)
    assert_that(cg_star >= cg)
    
    p = p + annotate('text', x = 0, y = -30, label = paste('cg*_', k, ' = ', cg_star, '\ncg_', k, ' = ', cg), hjust = 0, vjust = 1, size = 3)
  }
  
  if (add_district_labels){
    p = p + geom_sf_text(data = bogota, aes(label = LocNombre), colour = "black", size = 2)
  }
  p
  
}
