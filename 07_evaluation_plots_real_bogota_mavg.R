###############
# Plots for evaluation and hot spot predictions MAVG
# Nil-Jana Akpinar, 2020
###############

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(sf)

source('utils/real_bogota_utils.R')

set.seed(0)

bogota_shp = 'metadata/bogota.shp'
bogota = st_read(bogota_shp, stringsAsFactors = FALSE)

very_beginning = Sys.time()


################
# Adjustable parameters
################

version = 1
output_folder = paste0('version_', version)

# Parameters from data
time_steps = 365 * 6
assert_that(time_steps/(365/2) == round(time_steps/(365/2)))
cell_length = 1

# Evaluation time steps
eval_days = 2001:2190


################
# Assert that all integrals pre-computed
################

for (ts in eval_days) {
  
  true_integrals_path_ts = paste0('data/real_bogota_data_v', version, '_ts', time_steps, '/grid_integrals/grid_integrals_ts', ts, '.csv')
  full_estim_preds_path_ts = paste0('output/', output_folder, '/mavg_predictions/thinning_FALSE/grid_pred_ts', ts, '.csv')
  thinned_estim_preds_path_ts = paste0('output/', output_folder, '/mavg_predictions/thinning_TRUE/grid_pred_ts', ts, '.csv')
  
  if (!file.exists(true_integrals_path_ts)) {
    stop(paste0('True integrals for t = ', ts, ' not precomputed.'))
  }
  if (!file.exists(full_estim_preds_path_ts)) {
    stop(paste0('Full estimated integrals for t = ', ts, ' not precomputed.'))
  }
  if (!file.exists(thinned_estim_preds_path_ts)) {
    stop(paste0('Thinned estimated integrals for t = ', ts, ' not precomputed.'))
  }
  
}

###############
# Hotspot plots
###############

dir.create(paste0('output/', output_folder, '/mavg_hotspot_plots/'), showWarnings = FALSE)

q_full = list()
q_thinned = list()
i = 1

k = 50

full_thinned_estim_int_boxplot = data.frame(victim_true_rank = c(), 
                                            full_estim_rank = c(), 
                                            thinned_estim_rank = c(), 
                                            t = c())
full_thinned_estim_int_boxplot2 = data.frame(victim_int = c(), 
                                             victim_true_rank = c(),
                                             full_estim_rank = c(), 
                                             thinned_estim_rank = c(), 
                                             t = c())

for (ts in eval_days){
  
  # Read pre-computed integrals
  true_integrals_path_ts = paste0('data/real_bogota_data_v', version, '_ts', time_steps, '/grid_integrals/grid_integrals_ts', ts, '.csv')
  full_estim_preds_path_ts = paste0('output/', output_folder, '/mavg_predictions/thinning_FALSE/grid_pred_ts', ts, '.csv')
  thinned_estim_preds_path_ts = paste0('output/', output_folder, '/mavg_predictions/thinning_TRUE/grid_pred_ts', ts, '.csv')
  
  true_ints = read.csv(true_integrals_path_ts)
  full_estim_ints = read.csv(full_estim_preds_path_ts)
  thinned_estim_ints = read.csv(thinned_estim_preds_path_ts)
  
  # Make work with SEPP plotting function
  colnames(full_estim_ints)[colnames(full_estim_ints) == 'pred'] = 'int'
  colnames(thinned_estim_ints)[colnames(thinned_estim_ints) == 'pred'] = 'int'
  
  # Make hot spot plots over time
  q_full[[i]] = plot_hotspots(true_ints = true_ints,
                         estim_ints = full_estim_ints,
                         title = paste0('t = ', ts),
                         alpha = 0.4,
                         k = k,
                         add_cg = TRUE,
                         add_district_labels = FALSE)
  
  q_thinned[[i]] = plot_hotspots(true_ints = true_ints,
                              estim_ints = thinned_estim_ints,
                              title = paste0('t = ', ts),
                              alpha = 0.4,
                              k = k,
                              add_cg = TRUE,
                              add_district_labels = FALSE)
  
  # Look at how k cut off
  true_ints = true_ints[order(-true_ints$victim_int),]
  rownames(true_ints) = NULL
  true_ints$rank = 1:nrow(true_ints)
  q_true = ggplot(data = true_ints) + 
    geom_point(aes(x = rank, y = victim_int)) + 
    theme_bogota() + 
    geom_vline(xintercept = k) + 
    ylab('Hotspot integrals') + 
    xlab('Rank') + 
    ggtitle(paste0('True victimization integrals t = ', ts))

  full_estim_ints = full_estim_ints[order(-full_estim_ints$int),]
  rownames(full_estim_ints) = NULL
  full_estim_ints$rank = 1:nrow(full_estim_ints)
  q_full_estim = ggplot(data = full_estim_ints) + 
    geom_point(aes(x = rank, y = int)) + 
    theme_bogota() + 
    geom_vline(xintercept = k) + 
    ylab('Hotspot integrals') + 
    xlab('Rank') + 
    ggtitle(paste0('Predicted victimization integrals (no thinning) t = ', ts))

  thinned_estim_ints = thinned_estim_ints[order(-thinned_estim_ints$int),]
  rownames(thinned_estim_ints) = NULL
  thinned_estim_ints$rank = 1:nrow(thinned_estim_ints)
  q_thinned_estim = ggplot(data = thinned_estim_ints) + 
    geom_point(aes(x = rank, y = int)) + 
    theme_bogota() + 
    geom_vline(xintercept = k) + 
    ylab('Hotspot integrals') + 
    xlab('Rank') + 
    ggtitle(paste0('Predicted reporting integrals (thinning) t = ', ts))
  
  full_thinned_estim_int = full_estim_ints[,c('x', 'y', 'district', 'int', 'rank')]
  colnames(full_thinned_estim_int) = c('x', 'y', 'district', 'full_estim_int', 'full_estim_rank')
  full_thinned_estim_int = merge(x = full_thinned_estim_int, y = thinned_estim_ints[,c('x', 'y', 'district', 'int', 'rank')], all.x = TRUE)
  colnames(full_thinned_estim_int)[colnames(full_thinned_estim_int) == 'int'] = 'thinned_estim_int'
  colnames(full_thinned_estim_int)[colnames(full_thinned_estim_int) == 'rank'] = 'thinned_estim_rank'
  full_thinned_estim_int = merge(x = full_thinned_estim_int, y = true_ints[,c('x', 'y', 'district', 'victim_int', 'rank')])
  colnames(full_thinned_estim_int)[colnames(full_thinned_estim_int) == 'rank'] = 'victim_true_rank'
  
  q1 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = thinned_estim_rank, y = full_estim_rank, color = district)) + 
    geom_abline(intercept = 0, slope = 1)
  
  q2 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = victim_true_rank, y = full_estim_rank, color = district)) + 
    geom_abline(intercept = 0, slope = 1)
  
  q3 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = victim_true_rank, y = thinned_estim_rank, color = district)) + 
    geom_abline(intercept = 0, slope = 1)
  
  q4 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = full_estim_rank, y = victim_int, color = district))
  
  q5 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = thinned_estim_rank, y = victim_int, color = district))
  
  # Don't create hundreds of plots
  if (i <= 10){
    q = arrangeGrob(grobs = list(q_true, q_full_estim, q_thinned_estim, q1, q2, q3, q4, q5), nrow = 4, ncol = 3, top = paste0('True and predicted hotspots with k = ', k))
    pdf(paste0('output/', output_folder, '/mavg_hotspot_plots/hotspots_k', k, '_t', ts, '.pdf'), height = 20, width = 20)
    grid.draw(q)
    dev.off()
  }
  
  # Add to data frame for boxplots
  full_thinned_estim_int$t = ts
  full_thinned_estim_int_boxplot = rbind(full_thinned_estim_int_boxplot, full_thinned_estim_int[,c('victim_true_rank', 'full_estim_rank', 'thinned_estim_rank', 't')])
  full_thinned_estim_int_boxplot2 = rbind(full_thinned_estim_int_boxplot2, full_thinned_estim_int[,c('victim_int', 'victim_true_rank', 'full_estim_rank', 'thinned_estim_rank', 't')])
  
  i = i + 1
}

# Make boxplot
df_long2 = melt(full_thinned_estim_int_boxplot, id.vars = c('victim_true_rank', 't'))
p = ggplot(data = df_long2[df_long2$victim_true_rank <= 100,], aes(x = factor(victim_true_rank), y = value)) + 
  geom_boxplot() + 
  facet_wrap(~variable, ncol = 1) +
  theme_bogota() + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle(paste0('True hotspot ranks vs predicted hotspot ranks for top 100 true hotspots and \nt = ', min(eval_days), ',...,', max(eval_days)))

ggsave(filename = paste0('output/', output_folder, '/mavg_hotspot_plots/rank_boxplot_t', min(eval_days), '_', max(eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")

df_long2 = melt(full_thinned_estim_int_boxplot2, id.vars = c('victim_int', 't'))
p = ggplot(data = df_long2[(df_long2$value <= 100),], aes(x = factor(value), y = victim_int)) + 
  geom_boxplot() + 
  facet_wrap(~variable, ncol = 1) +
  theme_bogota() + 
  ggtitle(paste0('True crime in top 100 predicted hotspots, \nt = ', min(eval_days), ',...,', max(eval_days)))

ggsave(filename = paste0('output/', output_folder, '/mavg_hotspot_plots/int_boxplot_t', min(eval_days), '_', max(eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")


df_long2 = melt(full_thinned_estim_int_boxplot2, id.vars = c('victim_int', 't'))
p = ggplot(data = df_long2, aes(x = value, y = victim_int, color = variable)) + 
  stat_summary(geom = "line", 
               fun.y = mean) + 
  xlab('top {true, predicted with full data, predicted with reporting data} hotspots') + 
  ylab('True intensity in hotspot (i.e. crime found)') + 
  ggtitle(paste0('True crime in top 100 {true, predicted with full data, predicted with reporting data} hotspots \n(average over t = ', min(eval_days), ',...,', max(eval_days), ')')) + 
  xlim(c(0,100))

ggsave(filename = paste0('output/', output_folder, '/mavg_hotspot_plots/efficiency_gap_t', min(eval_days), '_', max(eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")


# Layout grid
idx_list = split(1:length(q_full), ceiling(seq_along(1:length(q_full))/12))
n_row = 2
for (idx in idx_list){
  
  idx = as.vector(unlist(idx))
  
  p_full = arrangeGrob(grobs = q_full[idx], nrow = n_row, top = paste0('Hotspots in evaluation period using victimization data (no thinning)'))
  p_thinned = arrangeGrob(grobs = q_thinned[idx], nrow = n_row, top = paste0('Hotspots in evaluation period using reportance data (thinning)'))
  
  # Write to pdf
  pdf(paste0('output/', output_folder, '/mavg_hotspot_plots/hotspots_full_k', k, '_tmin', min(eval_days) + min(idx) - 1, '_tmax', min(eval_days) + max(idx) - 1, '.pdf'), height = 15, width = 20)
  grid.draw(p_full)
  dev.off()
  
  pdf(paste0('output/', output_folder, '/mavg_hotspot_plots/hotspots_thinned_k', k, '_tmin', min(eval_days) + min(idx) - 1, '_tmax', min(eval_days) + max(idx) - 1, '.pdf'), height = 15, width = 20)
  grid.draw(p_thinned)
  dev.off()
}

###############

very_end = Sys.time()
time_delta = very_end - very_beginning
cat(paste0('\n\nRun time total: ', time_delta, ' ', units(time_delta), '\n'))

