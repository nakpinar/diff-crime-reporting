###############
# Plots for evaluation, sanity checks and hot spot predictions SEPP
# Nil-Jana Akpinar, 2020
###############

library(cubature)
library(spatstat)
library(plyr)
library(svMisc)
library(gridExtra)
library(reshape2)
library(fastDummies)
library(sf)

set.seed(0)

source('utils/real_bogota_utils.R')

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

# Parameters for evaluation
n_eval_days = 189
eval_days = (training_time_cutoff + 1):(training_time_cutoff + n_eval_days)
assert_that(max(eval_days) < time_steps)
assert_that(min(eval_days) > training_time_cutoff)


################
# Assert that all integrals pre-computed
################

for (ts in eval_days) {
  
  true_integrals_path_ts = paste0(true_integrals_path, 'grid_integrals_ts', ts, '.csv')
  full_estim_integrals_path_ts = paste0(full_estim_integrals_path, 'grid_integrals_ts', ts, '.csv')
  thinned_estim_integrals_path_ts = paste0(thinned_estim_integrals_path, 'grid_integrals_ts', ts, '.csv')
  
  if (!file.exists(true_integrals_path_ts)) {
    stop(paste0('True integrals for t = ', ts, ' not precomputed.'))
  }
  if (!file.exists(full_estim_integrals_path_ts)) {
    stop(paste0('Full estimated integrals for t = ', ts, ' not precomputed.'))
  }
  if (!file.exists(thinned_estim_integrals_path_ts)) {
    stop(paste0('Thinned estimated integrals for t = ', ts, ' not precomputed.'))
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
# Map plots vicitmization and reporting rates
###############

meta2 = meta
colnames(meta2)[colnames(meta2) == 'District'] = 'LocNombre'
bogota = merge(x = bogota, y = meta2, by = 'LocNombre', all.x = TRUE)

# Population density per km^2 (computed with grid cells)
bogota_mask = get_bogota_mask(cell_length = 1)
cells_per_district = data.frame(table(bogota_mask$district))
colnames(cells_per_district) = c('LocNombre', 'cells_length_1')
bogota = merge(x = bogota, y = cells_per_district, by = 'LocNombre', all.x = TRUE)
bogota$Population_density = bogota$Population / bogota$cells_length_1

# Absolute population in Bogota
pop1 = ggplot() + 
  geom_sf(data = bogota, aes(fill = Population)) + 
  theme_bogota() + 
  ggtitle('Absolute population') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

pop2 = ggplot() + 
  geom_sf(data = bogota, aes(fill = Population_density)) + 
  theme_bogota() + 
  ggtitle('Population density per km^2') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Number of victims
bogota$n_crimes_unscaled = bogota$Population * bogota$Victimization
victim2 = ggplot() + 
  geom_sf(data = bogota, aes(fill = n_crimes_unscaled)) + 
  theme_bogota() + 
  ggtitle('Absolute victimization') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Victimization rate
victim1 = ggplot() + 
  geom_sf(data = bogota, aes(fill = Victimization)) + 
  theme_bogota() + 
  ggtitle('Victimization rate (by population)') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Number of victims scaled by km^2
bogota$victim_area = bogota$n_crimes_unscaled / bogota$cells_length_1
victim3 = ggplot() + 
  geom_sf(data = bogota, aes(fill = victim_area)) + 
  theme_bogota() + 
  ggtitle('Absolute victimization per km^2') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Reporting rate
rep1 = ggplot() + 
  geom_sf(data = bogota, aes(fill = Percent_reported)) + 
  theme_bogota() + 
  ggtitle('Rate of crimes reported') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Crimes reported
bogota$abs_reported = bogota$n_crimes_unscaled * bogota$Percent_reported
rep2 = ggplot() + 
  geom_sf(data = bogota, aes(fill = abs_reported)) + 
  theme_bogota() + 
  ggtitle('Absolute crimes reported') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Crimes reported per area
bogota$abs_reported_area = bogota$abs_reported / bogota$cells_length_1
rep4 = ggplot() + 
  geom_sf(data = bogota, aes(fill = abs_reported_area)) + 
  theme_bogota() + 
  ggtitle('Absolute crimes reported per km^2') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

# Crimes reported per person
bogota$abs_reported_pop = bogota$abs_reported / bogota$Population
rep5 = ggplot() + 
  geom_sf(data = bogota, aes(fill = abs_reported_pop)) + 
  theme_bogota() + 
  ggtitle('Reported crimes per person') + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")


p_grid = arrangeGrob(grobs = list(pop1, pop2, victim1, victim2, victim3, rep1, rep2, rep5),
                     nrow = 3,
                     top = paste0('Bogota district crime reporting and victimization'))

# Write to file
pdf(paste0('mapplots/victim_rep_meta_data.pdf'), height = 20, width = 20)
grid.draw(p_grid)
dev.off()


###############
# Sanity checks: data and true intensity function vs survey meta data
###############

# Choose days for data check (should also span training data)
sanity_eval_days = 500:2190

# Create sanity checks output folder
sanity_path = paste0('data/real_bogota_data_v', version, '_ts', time_steps, '/sanity_check_plots/')
dir.create(sanity_path, showWarnings = FALSE)


## Map plots of crime implied by integrated intensities
# Victimization in meta, data and according to integrated intensities
true_integrals_path_ts = paste0(true_integrals_path, 'grid_integrals_ts', sanity_eval_days[1], '.csv')
true_int_mean = read.csv(true_integrals_path_ts)
colnames(true_int_mean)[colnames(true_int_mean) == 'victim_int'] = 'mean_int'
for (ts in sanity_eval_days[2:length(sanity_eval_days)]){
  true_integrals_path_ts = paste0(true_integrals_path, 'grid_integrals_ts', ts, '.csv')
  true_ints = read.csv(true_integrals_path_ts)[, c('x', 'y', 'district', 'victim_int')]
  true_int_mean = merge(x = true_int_mean, y = true_ints, all.x = TRUE)
  true_int_mean$mean_int = true_int_mean$mean_int + true_int_mean$victim_int
  true_int_mean = true_int_mean[,c('x', 'y', 'district', 'mean_int')]
}
true_int_mean$mean_int = true_int_mean$mean_int / length(sanity_eval_days)
true_int_mean$mean_int_semester = true_int_mean$mean_int * (365 / 2)

# Per cell 
p1 = ggplot() + 
  geom_bogota() + 
  theme_bogota() + 
  geom_raster(data = true_int_mean, aes(x = x, y = y, fill = mean_int), hjust = 0, vjust = 0) + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", name = 'crimes') +
  ggtitle(paste0('Daily crime as implied by true intensity integrals\n(average over t = ', min(sanity_eval_days), ',...,', max(sanity_eval_days), ')')) + 
  geom_bogota()

# Per district
true_int_mean_district = aggregate(true_int_mean$mean_int, by = list(Category = true_int_mean$district), FUN=sum)
colnames(true_int_mean_district) = c('LocNombre', 'n_crimes_integrals')
bogota = merge(x = bogota, y = true_int_mean_district, by = 'LocNombre', all.x = TRUE)

p2 = ggplot() + 
  geom_sf(data = bogota, aes(fill = n_crimes_integrals)) + 
  theme_bw() + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", name = 'crime') + 
  ggtitle(paste0('Daily crime as implied by true intensity integrals\n(average over t = ', min(sanity_eval_days), ',...,', max(sanity_eval_days), ')'))

p = arrangeGrob(p1, p2, nrow = 1)
ggsave(filename = paste0(sanity_path, 'daily_crime_true_integrated_intentisities_t', min(sanity_eval_days), '_', max(sanity_eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")


## Box plots of crime implied by true integrated intensities & crime implied by sampled data & desired crime by survey meta data 
# Get day-by-day sum of integrated intensities in districts
df = data.frame(district = c(), int_crime_t = c(), t = c())
for (ts in sanity_eval_days){
  
  true_integrals_path_ts = paste0(true_integrals_path, 'grid_integrals_ts', ts, '.csv')
  true_int = read.csv(true_integrals_path_ts)
  
  dd = aggregate(true_int$victim_int, by = list(Category = true_int$district), FUN=sum)
  colnames(dd) = c('district', 'int_crime_t')
  dd$t = ts
  
  df = rbind(df, dd)
}

# Get day-by-day sum of sampled crimes in districts
data_path = paste0('data/real_bogota_victim_data_v', version, '_ts', time_steps, '.csv')
victim_data = read.csv(data_path)
plot_data = victim_data[victim_data$t %in% sanity_eval_days,]
plot_data$count = 1
dd = aggregate(plot_data$count ~ plot_data$t + plot_data$district, FUN = sum)
colnames(dd) = c('t', 'district', 'data_crime_t')
dd = dummy_rows(dd, select_columns = c('t', 'district'), dummy_value = 0) # Insert missing 0 count rows

# Merge and plot
df = merge(x = df, y = dd, all.x = TRUE)
df_long = melt(df, id.vars=c('district', 't'))
rel_meta = meta[c('District', 'population_scaled', 'Victimization')]
colnames(rel_meta)[1] = 'district'
rel_meta$n_crimes_daily = rel_meta$population_scaled * rel_meta$Victimization / (365/2)
df_long = merge(x = df_long, y = rel_meta , by = c('district'), all.x = TRUE)
df_long$variable = revalue(df_long$variable, c("int_crime_t" = "integral", "data_crime_t" = "data")) # Rename factor values

p = ggplot(data = df_long, aes(x = variable, y = value)) + 
  geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 18, size = 2, color = "blue", fill = "blue") + 
  geom_hline(aes(yintercept = n_crimes_daily, linetype = 'CCB victimization \nsurvey crime'), color = 'red') +
  facet_wrap(~district, scales = 'free', ncol = 4) + 
  #ggtitle(paste0('Daily crime as implied by true intensity integrals (integral) and sampled victim data (data) \n(t = ', min(sanity_eval_days), ',...,', max(sanity_eval_days), ')')) + 
  theme_bw() + 
  xlab('') + 
  theme(legend.key = element_blank(), strip.background = element_rect(colour = "black", fill = "white")) + 
  ylab('Daily crime rate') + 
  scale_linetype_manual(name = '', values = c(2)) + 
  theme(legend.position = c(0.9,0.1))

p

ggsave(filename = paste0(sanity_path, 'boxplot_t', min(sanity_eval_days), '_', max(sanity_eval_days), '.png'),
       plot = p, width = 12, height = 14.5, dpi = 300, units = "in")


###############
# Hotspot plots
###############

dir.create(paste0('output/', output_folder, '/hotspot_plots/'), showWarnings = FALSE)

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
  true_ints = read.csv(paste0(true_integrals_path, 'grid_integrals_ts', ts, '.csv'))
  full_estim_ints = read.csv(paste0(full_estim_integrals_path, 'grid_integrals_ts', ts, '.csv'))
  thinned_estim_ints = read.csv(paste0(thinned_estim_integrals_path, 'grid_integrals_ts', ts, '.csv'))
  
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
  
  # Look at cut offs for k
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
  
  q4 = ggplot(data = full_thinned_estim_int[full_thinned_estim_int$thinned_estim_rank <= 100,]) + 
    geom_point(aes(x = thinned_estim_rank, y = full_estim_rank, color = district)) + 
    geom_abline(intercept = 0, slope = 1)
  
  q5 = ggplot(data = full_thinned_estim_int[full_thinned_estim_int$victim_true_rank <= 100,]) + 
    geom_point(aes(x = victim_true_rank, y = full_estim_rank, color = district)) + 
    geom_abline(intercept = 0, slope = 1)
  
  q6 = ggplot(data = full_thinned_estim_int[full_thinned_estim_int$victim_true_rank <= 100,]) + 
    geom_point(aes(x = victim_true_rank, y = thinned_estim_rank, color = district)) + 
    geom_abline(intercept = 0, slope = 1)
  
  q7 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = full_estim_rank, y = victim_int, color = district))
  
  q8 = ggplot(data = full_thinned_estim_int) + 
    geom_point(aes(x = thinned_estim_rank, y = victim_int, color = district))
  
  # Don't create hundreds of plots
  if (i <= 10){
    q = arrangeGrob(grobs = list(q_true, q_full_estim, q_thinned_estim, q1, q2, q3, q4, q5, q6, q7, q8), nrow = 4, ncol = 3, top = paste0('True and predicted hotspots with k = ', k))
    pdf(paste0('output/', output_folder, '/hotspot_plots/', 'hotspots_k', k, '_t', ts, '.pdf'), height = 20, width = 20)
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

ggsave(filename = paste0('output/', output_folder, '/hotspot_plots/rank_boxplot_t', min(eval_days), '_', max(eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")

df_long2 = melt(full_thinned_estim_int_boxplot2, id.vars = c('victim_int', 't'))
p = ggplot(data = df_long2[(df_long2$value <= 100),], aes(x = factor(value), y = victim_int)) + 
  geom_boxplot() + 
  facet_wrap(~variable, ncol = 1) +
  theme_bogota() + 
  ggtitle(paste0('True crime in top 100 predicted hotspots, \nt = ', min(eval_days), ',...,', max(eval_days)))

ggsave(filename = paste0('output/', output_folder, '/hotspot_plots/int_boxplot_t', min(eval_days), '_', max(eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")


df_long2 = melt(full_thinned_estim_int_boxplot2, id.vars = c('victim_int', 't'))
p = ggplot(data = df_long2, aes(x = value, y = victim_int, color = variable)) + 
  stat_summary(geom = "line", 
               fun.y = mean) + 
  xlab('top {true, predicted with full data, predicted with reporting data} hotspots') + 
  ylab('True intensity in hotspot (i.e. crime found)') + 
  ggtitle(paste0('True crime in top 100 {true, predicted with full data, predicted with reporting data} hotspots \n(average over t = ', min(eval_days), ',...,', max(eval_days), ')')) + 
  xlim(c(0,100))

ggsave(filename = paste0('output/', output_folder, '/hotspot_plots/efficiency_gap_t', min(eval_days), '_', max(eval_days), '.png'),
       plot = p, width = 10, height = 8, dpi = 300, units = "in")

# Layout grid (each 12 on new grid)
idx_list = split(1:length(q_full), ceiling(seq_along(1:length(q_full))/12))
n_row = 2

for (idx in idx_list){
  
  idx = as.vector(unlist(idx))
  
  p_full = arrangeGrob(grobs = q_full[idx], nrow = n_row, top = paste0('Hotspots in evaluation period using victimization data (no thinning)'))
  p_thinned = arrangeGrob(grobs = q_thinned[idx], nrow = n_row, top = paste0('Hotspots in evaluation period using reportance data (thinning)'))
  
  # Write to pdf
  pdf(paste0('output/', output_folder, '/hotspot_plots/', 'hotspots_full_k', k, '_tmin', min(eval_days) + min(idx) - 1, '_tmax', min(eval_days) + max(idx) - 1, '.pdf'), height = 15, width = 20)
  grid.draw(p_full)
  dev.off()
  
  pdf(paste0('output/', output_folder, '/hotspot_plots/', 'hotspots_thinned_k', k, '_tmin', min(eval_days) + min(idx) - 1, '_tmax', min(eval_days) + max(idx) - 1, '.pdf'), height = 15, width = 20)
  grid.draw(p_thinned)
  dev.off()
}

###############

very_end = Sys.time()
time_delta = very_end - very_beginning
cat(paste0('\n\nRun time total: ', time_delta, ' ', units(time_delta), '\n'))

