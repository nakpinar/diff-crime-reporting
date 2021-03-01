################
# District fairness plots over all models; Plots are created aggregate over simulation runs (versions).
# Nil-Jana Akpinar, 2020
################

library(cubature)
library(spatstat)
library(plyr)
library(svMisc)
library(gridExtra)
library(reshape2)
library(fastDummies)
library(sf)
library(assertthat)

set.seed(0)

source('utils/real_bogota_utils.R')

bogota_shp = 'metadata/bogota.shp'
bogota = st_read(bogota_shp, stringsAsFactors = FALSE)

very_beginning = Sys.time()


###############
# Adjustable parameters
###############

# All versions to aggregate over for plot
versions = c(1)

# Number of hot spots
k = 50


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

colnames(meta)[colnames(meta) == 'District'] = 'district'

###############
# Measure district fairness function
###############

dir.create(paste0('output/version_', min(versions), '_to_', max(versions), '/'), showWarnings = FALSE)

# Returns data frame with relevant counts and cg for given time step
get_fairness_districts = function(true_ints, full_estim_ints, thinned_estim_ints, full_estim_preds, thinned_estim_preds, ts, k){
  
  # Add mask for hot spots
  true_ints = true_ints[order(-true_ints$victim_int),]
  rownames(true_ints) = NULL
  true_ints$int_true_topk = 0
  true_ints$int_true_topk[1:k] = true_ints$victim_int[1:k]
  true_ints$true_topk = 0
  true_ints$true_topk[1:k] = 1
  
  full_estim_ints = full_estim_ints[order(-full_estim_ints$int),] 
  rownames(full_estim_ints) = NULL
  full_estim_ints$int_full_topk = 0
  full_estim_ints$int_full_topk[1:k] = full_estim_ints$int[1:k]
  full_estim_ints$full_topk = 0
  full_estim_ints$full_topk[1:k] = 1
  
  thinned_estim_ints = thinned_estim_ints[order(-thinned_estim_ints$int),] 
  rownames(thinned_estim_ints) = NULL
  thinned_estim_ints$int_thinned_topk = 0
  thinned_estim_ints$int_thinned_topk[1:k] = thinned_estim_ints$int[1:k]
  thinned_estim_ints$thinned_topk = 0
  thinned_estim_ints$thinned_topk[1:k] = 1
  
  full_estim_preds = full_estim_preds[order(-full_estim_preds$pred),] 
  rownames(full_estim_preds) = NULL
  full_estim_preds$pred_full_topk = 0
  full_estim_preds$pred_full_topk[1:k] = full_estim_ints$int[1:k]
  full_estim_preds$mavg_full_topk = 0
  full_estim_preds$mavg_full_topk[1:k] = 1
  colnames(full_estim_preds)[colnames(full_estim_preds) == 'pred'] = 'full_pred'
  
  thinned_estim_preds = thinned_estim_preds[order(-thinned_estim_preds$pred),] 
  rownames(thinned_estim_preds) = NULL
  thinned_estim_preds$pred_thinned_topk = 0
  thinned_estim_preds$pred_thinned_topk[1:k] = thinned_estim_ints$int[1:k]
  thinned_estim_preds$mavg_thinned_topk = 0
  thinned_estim_preds$mavg_thinned_topk[1:k] = 1
  colnames(thinned_estim_preds)[colnames(thinned_estim_preds) == 'pred'] = 'thinned_pred'
  
  # Rescale thinned estimated intensities by rate of underreporting
  thinned_estim_ints_rescaled = thinned_estim_ints[,c('x', 'y', 'district', 'cell_id', 'int')]
  thinned_estim_ints_rescaled = merge(x = thinned_estim_ints_rescaled, y = meta[,c('district', 'Percent_reported')], all.x = TRUE)
  thinned_estim_ints_rescaled$rescaled_int = thinned_estim_ints_rescaled$int / thinned_estim_ints_rescaled$Percent_reported
  thinned_estim_ints_rescaled = thinned_estim_ints_rescaled[order(-thinned_estim_ints_rescaled$rescaled_int),] 
  rownames(thinned_estim_ints_rescaled) = NULL
  thinned_estim_ints_rescaled$rescaled_int_thinned_topk = 0
  thinned_estim_ints_rescaled$rescaled_int_thinned_topk[1:k] = thinned_estim_ints_rescaled$rescaled_int[1:k]
  thinned_estim_ints_rescaled$rescaled_thinned_topk = 0
  thinned_estim_ints_rescaled$rescaled_thinned_topk[1:k] = 1
  
  thinned_estim_preds_rescaled = thinned_estim_preds[,c('x', 'y', 'district', 'cell_id', 'thinned_pred')]
  thinned_estim_preds_rescaled = merge(x = thinned_estim_preds_rescaled, y = meta[,c('district', 'Percent_reported')], all.x = TRUE)
  thinned_estim_preds_rescaled$rescaled_preds = thinned_estim_preds_rescaled$thinned_pred / thinned_estim_preds_rescaled$Percent_reported
  thinned_estim_preds_rescaled = thinned_estim_preds_rescaled[order(-thinned_estim_preds_rescaled$rescaled_preds),] 
  rownames(thinned_estim_preds_rescaled) = NULL
  thinned_estim_preds_rescaled$rescaled_preds_thinned_topk = 0
  thinned_estim_preds_rescaled$rescaled_preds_thinned_topk[1:k] = thinned_estim_preds_rescaled$rescaled_preds[1:k]
  thinned_estim_preds_rescaled$rescaled_mavg_thinned_topk = 0
  thinned_estim_preds_rescaled$rescaled_mavg_thinned_topk[1:k] = 1
  
  # Get counts, cgs, cg/count
  all = merge(x = true_ints, y = full_estim_ints, by = c('x', 'y', 'district'), all.x = TRUE)
  all = merge(x = all, y = thinned_estim_ints, by = c('x', 'y', 'district'), all.x = TRUE)
  all = merge(x = all, y = full_estim_preds, by = c('x', 'y', 'district', 'cell_id'), all.x = TRUE)
  all = merge(x = all, y = thinned_estim_preds, by = c('x', 'y', 'district', 'cell_id'), all.x = TRUE)
  all = merge(x = all, y = thinned_estim_ints_rescaled, by = c('x', 'y', 'district', 'cell_id'), all.x = TRUE)
  all = merge(x = all, y = thinned_estim_preds_rescaled, by = c('x', 'y', 'district', 'cell_id'), all.x = TRUE)
  
  counts = aggregate(list(true_topk = all$true_topk,
                          full_topk = all$full_topk,
                          thinned_topk = all$thinned_topk,
                          mavg_full_topk = all$mavg_full_topk,
                          mavg_thinned_topk = all$mavg_thinned_topk,
                          rescaled_thinned_topk = all$rescaled_thinned_topk,
                          rescaled_mavg_thinned_topk = all$rescaled_mavg_thinned_topk),
                     by = list(district = all$district), FUN = sum)
  all$true_int_full_topk = all$full_topk * all$victim_int
  all$true_int_thinned_topk = all$thinned_topk * all$victim_int
  all$true_int_mavg_full_topk = all$mavg_full_topk * all$victim_int
  all$true_int_mavg_thinned_topk = all$mavg_thinned_topk * all$victim_int
  all$true_int_rescaled_thinned_topk = all$rescaled_thinned_topk * all$victim_int
  all$true_int_rescaled_mavg_thinned_topk = all$rescaled_mavg_thinned_topk * all$victim_int
  cgs = aggregate(list(cg_star = all$int_true_topk,
                       cg_full = all$true_int_full_topk,
                       cg_thinned = all$true_int_thinned_topk, 
                       cg_full_mavg = all$true_int_mavg_full_topk, 
                       cg_thinned_mavg = all$true_int_mavg_thinned_topk, 
                       cg_rescaled_thinned = all$true_int_rescaled_thinned_topk, 
                       cg_rescaled_thinned_mavg = all$true_int_rescaled_mavg_thinned_topk), 
                  by = list(district = all$district), FUN = sum)
  
  rel_counts = counts
  denominator = rel_counts$true_topk
  rel_counts[,2:ncol(rel_counts)] = rel_counts[,2:ncol(rel_counts)] / denominator
  # Replace Nan by 1000 to filter out later
  rel_counts[is.na(rel_counts)] = 1000
  # Replace Inf by original value * -1 to filter out later
  rel_counts = do.call(data.frame,lapply(rel_counts, function(x) replace(x, is.infinite(x),NA)))
  rel_counts[is.na(rel_counts)] = as.numeric(counts[is.na(rel_counts)]) * -1
  
  diff_counts = counts
  diff_counts[,2:ncol(diff_counts)] = diff_counts[,2:ncol(diff_counts)] - diff_counts$true_topk
  
  # Average cg, i.e. cg/count
  avg = merge(x = counts, y = cgs, by = c('district'), all.x = TRUE)
  avg$avg_cg_star = avg$cg_star / avg$true_topk
  avg$avg_cg_star[is.na(avg$avg_cg_star)] = 0
  avg$avg_cg_full = avg$cg_full / avg$full_topk
  avg$avg_cg_full[is.na(avg$avg_cg_full)] = 0
  avg$avg_cg_thinned = avg$cg_thinned / avg$thinned_topk
  avg$avg_cg_thinned[is.na(avg$avg_cg_thinned)] = 0
  avg$avg_cg_full_mavg = avg$cg_full_mavg / avg$mavg_full_topk
  avg$avg_cg_full_mavg[is.na(avg$avg_cg_full_mavg)] = 0
  avg$avg_cg_thinned_mavg = avg$cg_thinned_mavg / avg$mavg_thinned_topk
  avg$avg_cg_thinned_mavg[is.na(avg$avg_cg_thinned_mavg)] = 0
  avg$avg_cg_rescaled = avg$cg_rescaled_thinned / avg$rescaled_thinned_topk
  avg$avg_cg_rescaled[is.na(avg$avg_cg_rescaled)] = 0
  avg$avg_cg_rescaled_mavg = avg$cg_rescaled_thinned_mavg / avg$rescaled_mavg_thinned_topk
  avg$avg_cg_rescaled_mavg[is.na(avg$avg_cg_rescaled_mavg)] = 0
  avg = avg[, c('district', 'avg_cg_star', 'avg_cg_full', 'avg_cg_thinned', 'avg_cg_full_mavg', 'avg_cg_thinned_mavg', 'avg_cg_rescaled', 'avg_cg_rescaled_mavg')]
  
  avg2 = avg
  denominator = avg2$avg_cg_star
  denominator[denominator == 0] = 1
  avg2$avg_cg_full = avg2$avg_cg_full / denominator
  avg2$avg_cg_full[is.na(avg2$avg_cg_full)] = 0
  avg2$avg_cg_thinned = avg2$avg_cg_thinned / denominator
  avg2$avg_cg_thinned[is.na(avg2$avg_cg_thinned)] = 0
  avg2$avg_cg_full_mavg = avg2$avg_cg_full_mavg / denominator
  avg2$avg_cg_full_mavg[is.na(avg2$avg_cg_full_mavg)] = 0
  avg2$avg_cg_thinned_mavg = avg2$avg_cg_thinned_mavg / denominator
  avg2$avg_cg_thinned_mavg[is.na(avg2$avg_cg_thinned_mavg)] = 0
  avg2$avg_cg_star = avg2$avg_cg_star / denominator
  avg2$avg_cg_star[is.na(avg2$avg_cg_star)] = 0
  avg2$avg_cg_rescaled = avg2$avg_cg_rescaled / denominator
  avg2$avg_cg_rescaled[is.na(avg2$avg_cg_rescaled)] = 0
  avg2$avg_cg_rescaled_mavg = avg2$avg_cg_rescaled_mavg / denominator
  avg2$avg_cg_rescaled_mavg[is.na(avg2$avg_cg_rescaled_mavg)] = 0
  
  # Minimum cg* in any of the selected hotspots in the region
  all2 = all
  all2[all2 == 0] = 9999
  min_cg_star = aggregate(list(cg_star = all2$int_true_topk,
                               cg_full = all2$true_int_full_topk,
                               cg_thinned = all2$true_int_thinned_topk,
                               cg_full_mavg = all2$true_int_mavg_full_topk,
                               cg_thinned_mavg = all2$true_int_mavg_thinned_topk,
                               cg_rescaled_thinned = all2$true_int_rescaled_thinned_topk,
                               cg_rescaled_thinned_mavg = all2$true_int_rescaled_mavg_thinned_topk),
                          by = list(district = all2$district), FUN = min)
  
  # Sanity check
  assert_that(sum(cgs$cg_thinned) <= sum(cgs$cg_star))
  assert_that(sum(cgs$cg_full) <= sum(cgs$cg_star))
  assert_that(sum(cgs$cg_thinned_mavg) <= sum(cgs$cg_star))
  assert_that(sum(cgs$cg_full_mavg) <= sum(cgs$cg_star))
  
  #res = merge(x = counts, y = cgs, by = c('district'), all.x = TRUE)
  list(counts, cgs, avg, avg2, min_cg_star, rel_counts, diff_counts)
}


###############
# Iterate over versions to get fairness summary data
###############

all_version_rel_counts_long = data.frame(district = c(), t = c(), value = c(), version = c())
all_version_min_cg_long = data.frame(district = c(), t = c(), value = c(), version = c())
all_version_diff_counts_long = data.frame(district = c(), t = c(), value = c(), version = c())
all_version_counts_long = data.frame(district = c(), t = c(), value = c(), version = c())

for (version in versions){
  
cat('Start version:', version, '\n')

output_folder = paste0('version_', version)

# Parameters from data
time_steps = 365 * 6
assert_that(time_steps/(365/2) == round(time_steps/(365/2)))
background_sd = 15

# Path to integrated true victimization intensities
true_integrals_path = paste0('data/real_bogota_data_v', version, '_ts', time_steps, '/grid_integrals/')

# Time steps for training
ignore_first_t = 500 
training_time_cutoff = 4 * 500
cell_length = 1

# Path to integrated model intensities
full_estim_integrals_path = paste0('output/', output_folder, '/estim_grid_integrals/thinning_FALSE/')
thinned_estim_integrals_path = paste0('output/', output_folder, '/estim_grid_integrals/thinning_TRUE/')

# Path to MAVG predictions
full_estim_preds_path = paste0('output/', output_folder, '/mavg_predictions/thinning_FALSE/')
thinned_estim_preds_path = paste0('output/', output_folder, '/mavg_predictions/thinning_TRUE/')

# Parameters for evaluation
n_eval_days = 189
eval_days = (training_time_cutoff + 1):(training_time_cutoff + n_eval_days)
assert_that(max(eval_days) < time_steps)
assert_that(min(eval_days) > training_time_cutoff)

# Initiate data frames
counts_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
cgs_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
avg_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
avg2_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
min_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
min_long_all2 = data.frame(district = c(), t = c(), variable = c(), value = c())
rel_counts_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
diff_counts_long_all = data.frame(district = c(), t = c(), variable = c(), value = c())
  
# We want the number of predicted hotspots / number of true hotspots in each district for every time steps
for (ts in eval_days){
    
    # Read pre-computed integrals
    true_ints = read.csv(paste0(true_integrals_path, 'grid_integrals_ts', ts, '.csv'))
    full_estim_ints = read.csv(paste0(full_estim_integrals_path, 'grid_integrals_ts', ts, '.csv'))
    thinned_estim_ints = read.csv(paste0(thinned_estim_integrals_path, 'grid_integrals_ts', ts, '.csv'))
    full_estim_preds = read.csv(paste0(full_estim_preds_path, 'grid_pred_ts', ts, '.csv'))
    thinned_estim_preds = read.csv(paste0(thinned_estim_preds_path, 'grid_pred_ts', ts, '.csv'))
    
    # Get fairness counts and cgs
    res = get_fairness_districts(true_ints = true_ints,
                                 full_estim_ints = full_estim_ints,
                                 thinned_estim_ints = thinned_estim_ints,
                                 full_estim_preds = full_estim_preds,
                                 thinned_estim_preds = thinned_estim_preds,
                                 ts = ts,
                                 k = k)
    
    counts = res[[1]]
    cgs = res[[2]]
    avg = res[[3]]
    avg2 = res[[4]]
    min_cg_star = res[[5]]
    rel_counts = res[[6]]
    diff_counts = res[[7]]
    
    # Add to main data frame
    counts$t = ts
    counts_long = melt(counts, id.vars = c('district', 't'))
    counts_long_all = rbind(counts_long_all, counts_long)
    
    cgs$t = ts
    cgs_long = melt(cgs, id.vars = c('district', 't'))
    cgs_long_all = rbind(cgs_long_all, cgs_long)
    
    avg$t = ts
    avg_long = melt(avg, id.vars = c('district', 't'))
    avg_long_all = rbind(avg_long_all, avg_long)
    
    avg2$t = ts
    avg2_long = melt(avg2, id.vars = c('district', 't'))
    avg2_long_all = rbind(avg2_long_all, avg2_long)
    
    min_cg_star$t = ts
    min_long = melt(min_cg_star, id.vars = c('district', 't'))
    min_long_all = rbind(min_long_all, min_long)
    
    rel_counts$t = ts
    rel_counts_long = melt(rel_counts, id.vars = c('district', 't'))
    rel_counts_long_all = rbind(rel_counts_long_all, rel_counts_long)
    
    diff_counts$t = ts
    diff_counts_long = melt(diff_counts, id.vars = c('district', 't'))
    diff_counts_long_all = rbind(diff_counts_long_all, diff_counts_long)
    
    # Same, but we set 0 to cg_star_50 which is the best case for the model; ie instead column-wise min
    min_cg_star2 = min_cg_star
    for (i in 1:ncol(min_cg_star2)){
      if(grepl('cg', colnames(min_cg_star2)[i], fixed = TRUE)){
        min_cg_star2[,i][min_cg_star2[,i] == 0] = min(min_cg_star2[,i][min_cg_star2[,i] != 0])
      }
    }
    min_long2 = melt(min_cg_star2, id.vars = c('district', 't'))
    min_long_all2 = rbind(min_long_all2, min_long2)
    
}

# Keep track of min and rel_cg for all versions
min_long_all$version = version
all_version_min_cg_long = rbind(all_version_min_cg_long, min_long_all)

rel_counts_long_all$version = version
all_version_rel_counts_long = rbind(all_version_rel_counts_long, rel_counts_long_all)

diff_counts_long_all$version = version
all_version_diff_counts_long = rbind(all_version_diff_counts_long, diff_counts_long_all)

counts_long_all$version = version
all_version_counts_long = rbind(all_version_counts_long, counts_long_all)

} 

bogota_mask = get_bogota_mask(cell_length)
bogota_mask_help = bogota_mask[bogota_mask$district != 0,]

assert_that(nrow(all_version_rel_counts_long)/19/length(versions)/189 == 7)
assert_that(nrow(all_version_min_cg_long)/19/length(versions)/189 == 7)
assert_that(nrow(all_version_diff_counts_long)/19/length(versions)/189 == 7)
assert_that(nrow(all_version_counts_long)/19/length(versions)/189 == 7)


###############
# District fairness plots
###############

counts_wo_x_div_0 = all_version_rel_counts_long[all_version_rel_counts_long$value >= 0,]
counts_wo_x_div_0$value[counts_wo_x_div_0$value == 1000] = 1 # set 0/0 to 1

## Relative count
all_version_rel_counts_long$variable = revalue(all_version_rel_counts_long$variable, c('true_topk' = 'True', 'full_topk' = 'S1', 'thinned_topk' = 'S2', 'rescaled_thinned_topk' = 'S3', 'mavg_full_topk' = 'M1', 'mavg_thinned_topk' = 'M2', 'rescaled_mavg_thinned_topk' = 'M3'))
all_version_rel_counts_long$variable = factor(all_version_rel_counts_long$variable, levels = c('True', 'S1', 'M1', 'S3', 'M3', 'S2', 'M2'))
counts_wo_x_div_0 = all_version_rel_counts_long[all_version_rel_counts_long$value >= 0,]
counts_wo_x_div_0$value[counts_wo_x_div_0$value == 1000] = 1 # set 0/0 to 1
counts_wo_x_div_0 = counts_wo_x_div_0[counts_wo_x_div_0$variable != 'True th',]

p = ggplot(data = counts_wo_x_div_0[counts_wo_x_div_0$variable != 'True',], aes(x = factor(variable), y = value)) +
  facet_wrap(~district, scales = 'free', ncol = 4) + 
  theme_bw() +
  geom_boxplot(alpha = 1, outlier.size = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "blue", fill = "blue", shape = 18) +
  xlab('') + 
  ylab('(# Predicted hot spots in district) / (# True hot spots in district)') + 
  theme(legend.key = element_blank(), strip.background = element_rect(colour = "black", fill = "white")) + 
  ggtitle(paste0('(# Selected hotspots)/(# True hotspots) without x/0 but with 0/0 = 1, ', length(eval_days), ' evaluation days x ', length(versions), ' simulation runs, choosing top ', k, ' hotspots'))

ggsave(filename = paste0('output/version_', min(versions), '_to_', max(versions), '/rel_count.png'),
       plot = p, width = 12, height = 14.5, dpi = 300, units = "in")

## Difference in count
all_version_diff_counts_long$variable = revalue(all_version_diff_counts_long$variable, c('true_topk' = 'True', 'full_topk' = 'S1', 'thinned_topk' = 'S2', 'rescaled_thinned_topk' = 'S3', 'mavg_full_topk' = 'M1', 'mavg_thinned_topk' = 'M2', 'rescaled_mavg_thinned_topk' = 'M3'))
all_version_diff_counts_long$variable = factor(all_version_diff_counts_long$variable, levels = c('True', 'S1', 'M1', 'S3', 'M3', 'S2', 'M2'))

p = ggplot(data = all_version_diff_counts_long[all_version_diff_counts_long$variable != 'True',], aes(x = factor(variable), y = value)) +
  facet_wrap(~district, scales = 'free', ncol = 4) + 
  theme_bw() +
  geom_boxplot(alpha = 1, outlier.size = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "blue", fill = "blue", shape = 18) +
  xlab('') + 
  ylab('(# Predicted hot spots in district) - (# True hot spots in district)') + 
  theme(legend.key = element_blank(), strip.background = element_rect(colour = "black", fill = "white"))

ggsave(filename = paste0('output/version_', min(versions), '_to_', max(versions), '/diff_count.png'),
       plot = p, width = 12, height = 14.5, dpi = 300, units = "in")

## Relative count -- cases with no true hot spots
label1 = '(# Predicted hot spots) > 0'
label2 = '(# Predicted hot spots) = 0'

counts_x_div_0 =  all_version_rel_counts_long[all_version_rel_counts_long$value < 0,]
counts_x_div_0$case = label1
counts_0_div_0 =  all_version_rel_counts_long[all_version_rel_counts_long$value == 1000,]
counts_0_div_0$case = label2
counts_special = rbind(counts_0_div_0, counts_x_div_0)                                       
counts_special$count = 1

cases = aggregate(list(count = counts_special$count), by = list(district = counts_special$district,
                                                                variable = counts_special$variable,
                                                                case = counts_special$case),
                  FUN = sum)



cases$rel_count = cases$count / (length(versions)*189)

p = ggplot(data = cases[cases$variable != 'True',]) + 
  geom_bar(aes(x = factor(variable), y = rel_count, fill = case), alpha = 0.7, stat='identity') + 
  facet_wrap(~district, scales = 'free', ncol = 4) + 
  theme_bw() + 
  theme(legend.key = element_blank(), strip.background = element_rect(colour = "black", fill = "white")) + 
  ylab('Fraction of cases') + 
  xlab('') + 
  theme(legend.position = c(0.375,0.1)) +
  scale_fill_manual("", values = c( '(# Predicted hot spots) > 0' = "orange", '(# Predicted hot spots) = 0' = "blue"))

ggsave(filename = paste0('output/version_', min(versions), '_to_', max(versions), '/rel_count_no_true.png'),
       plot = p, width = 12, height = 14.5, dpi = 300, units = "in")

## Minimum CG
all_version_min_cg_long$variable = revalue(all_version_min_cg_long$variable, c('cg_star' = 'True', 'cg_full' = 'S1', 'cg_thinned' = 'S2', 'cg_rescaled_thinned' = 'S3', 'cg_full_mavg' = 'M1', 'cg_thinned_mavg' = 'M2', 'cg_rescaled_thinned_mavg' = 'M3'))
all_version_min_cg_long$variable = factor(all_version_min_cg_long$variable, levels = c('True', 'S1', 'M1', 'S3', 'M3', 'S2', 'M2'))
min_cgs = all_version_min_cg_long[all_version_min_cg_long$value != 9999,]

p = ggplot(data = min_cgs[min_cgs$variable != 'True',], aes(x = factor(variable), y = value)) +
  facet_wrap(~district, scales = 'free', ncol = 4) + 
  theme_bw() +
  geom_boxplot(alpha = 1, outlier.size = 0.1) + 
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "blue", fill = "blue", shape = 18) +
  xlab('') + 
  ylab('Mininum true crime rate in predicted hot spot') +
  theme(legend.key = element_blank(), strip.background = element_rect(colour = "black", fill = "white"))

ggsave(filename = paste0('output/version_', min(versions), '_to_', max(versions), '/min_cgs.png'),
       plot = p, width = 12, height = 14.5, dpi = 300, units = "in")