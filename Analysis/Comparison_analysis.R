rm(list = ls())
library(Seurat)
library(hdWGCNA)
library(diffuStats)
library(tidyverse)
library(igraph)
library(devtools)
# devtools::install_github("cysouw/qlcMatrix")
library(qlcMatrix)
library(ggpubr)
library(data.table)
gc()

load(file = './073124_data_to_anal.Rdata')
########################################################################################################################################

Thyroid_corr_1 <- MC_res %>%
  dplyr::filter(., Phenotype != 'hyper') %>%
  pivot_wider(id_cols = c('Module', 'gene', 'Cell_type'),names_from = Phenotype, values_from = Final.Heat) %>%
  group_by(Module, Cell_type) %>%
  dplyr::summarise(correlation.propagation = cor(Graves, hypo, use = 'pairwise.complete.obs', method = 'spearman')) %>%
  na.omit


Thyroid_corr_2 <- MC_res %>%
  dplyr::filter(., Phenotype != 'hypo') %>%
  pivot_wider(id_cols = c('Module', 'gene', 'Cell_type'),names_from = Phenotype, values_from = Final.Heat) %>%
  group_by(Module, Cell_type) %>%
  dplyr::summarise(correlation.propagation = cor(Graves, hyper, use = 'pairwise.complete.obs', method = 'spearman')) %>%
  na.omit



Thyroid_corr_TWAS_1 <- MC_res %>%
  dplyr::filter(., Phenotype != 'hyper') %>%
  dplyr::select(., -c(Final.Heat, delta)) %>%
  pivot_wider(id_cols = c('Module', 'gene', 'Cell_type'),names_from = Phenotype, values_from = TWAS.Z) %>%
  group_by(Module, Cell_type) %>%
  summarise(correlation.TWAS = cor(Graves, hypo, use = 'pairwise.complete.obs', method = 'spearman')) %>%
  na.omit

Thyroid_corr_TWAS_2 <- MC_res %>%
  dplyr::filter(., Phenotype != 'hypo') %>%
  dplyr::select(., -c(Final.Heat, delta)) %>%
  pivot_wider(id_cols = c('Module', 'gene', 'Cell_type'),names_from = Phenotype, values_from = TWAS.Z) %>%
  group_by(Module, Cell_type) %>%
  summarise(correlation.TWAS = cor(Graves, hyper, use = 'pairwise.complete.obs', method = 'spearman')) %>%
  na.omit


Thyroid_corr_delta_1 <- MC_res %>%
  dplyr::filter(., Phenotype != 'hyper') %>%
  dplyr::select(., -c(Final.Heat, TWAS.Z)) %>%
  pivot_wider(id_cols = c('Module', 'gene', 'Cell_type'), names_from = Phenotype, values_from = delta) %>%
  group_by(Module, Cell_type) %>%
  summarise(correlation.delta = cor(Graves, hypo, use = 'pairwise.complete.obs', method = 'spearman')) %>%
  na.omit


Thyroid_corr_delta_2 <- MC_res %>%
  dplyr::filter(., Phenotype != 'hypo') %>%
  dplyr::select(., -c(Final.Heat, TWAS.Z)) %>%
  pivot_wider(id_cols = c('Module', 'gene', 'Cell_type'), names_from = Phenotype, values_from = delta) %>%
  group_by(Module, Cell_type) %>%
  summarise(correlation.delta = cor(Graves, hyper, use = 'pairwise.complete.obs', method = 'spearman')) %>%
  na.omit

Thyroid_corr_1$Module %>% str

color_df_1 <- data.frame(Factor = unique(Thyroid_corr_1$Module),
                       Char = as.character(unique(Thyroid_corr_1$Module)))

color_df_2 <- data.frame(Factor = unique(Thyroid_corr_2$Module),
                         Char = as.character(unique(Thyroid_corr_2$Module)))


# Thyroid_corr$correlation.propagation %>% quantile(., c(0.1, 0.9)) %>% .[1]

my.quantile <- c(0.1, 0.9)

if(!dir.exists('./Grid_anal_073124')){
  dir.create('./Grid_anal_073124')
}


ggplot(Thyroid_corr_1, aes(x = Module, y = correlation.propagation)) +
  geom_col(aes(fill = Module), colour = 'grey70') +
  geom_hline(yintercept = quantile(Thyroid_corr_1$correlation.propagation, my.quantile, na.rm = T)[2], color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = quantile(Thyroid_corr_1$correlation.propagation, my.quantile, na.rm = T)[1], color = 'darkblue', linetype = 'dashed') +
  # labs(title = "Correlation between Hyperthyroidism and Hypothyroidism") +
  ylab("Spearman's rho") + 
  theme_bw() +
  scale_fill_manual(values = color_df_1$Char) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6)) +
  facet_wrap(~Cell_type, nrow = 2, scales = 'free_x')

ggsave('./Grid_anal_073124/Correlation_Final_heat_hypo_GD.png', width = 10, height = 3.5, dpi = 300, units = 'in')

ggplot(Thyroid_corr_2, aes(x = Module, y = correlation.propagation)) +
  geom_col(aes(fill = Module), colour = 'grey70') +
  geom_hline(yintercept = quantile(Thyroid_corr_2$correlation.propagation, my.quantile, na.rm = T)[2], color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = quantile(Thyroid_corr_2$correlation.propagation, my.quantile, na.rm = T)[1], color = 'darkblue', linetype = 'dashed') +
  # labs(title = "Correlation between Hyperthyroidism and Hypothyroidism") +
  ylab("Spearman's rho") + 
  theme_bw() +
  scale_fill_manual(values = color_df_2$Char)  +
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6)) +
  facet_wrap(~Cell_type, nrow = 2, scales = 'free_x')

ggsave('./Grid_anal_073124/Correlation_Final_heat_hyper_GD.png', width = 10, height = 3.5, dpi = 300, units = 'in')

ggplot(Thyroid_corr_TWAS_1, aes(x = Module, y = correlation.TWAS)) +
  geom_col(aes(fill = Module), colour = 'grey70') +
  geom_hline(yintercept = quantile(Thyroid_corr_TWAS_1$correlation.TWAS, my.quantile, na.rm = T)[2], color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = quantile(Thyroid_corr_TWAS_1$correlation.TWAS, my.quantile, na.rm = T)[1], color = 'darkblue', linetype = 'dashed') +
  # labs(title = "Correlation between Hyperthyroidism and Hypothyroidism") +
  ylab("Spearman's rho") + 
  theme_bw() +
  scale_fill_manual(values = color_df_1$Char) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6)) +
  facet_wrap(~Cell_type, nrow = 2, scales = 'free_x')

ggsave('./Grid_anal_073124/Correlation_TWAS_hypo_GD.png', width = 10, height = 3.5, dpi = 300, units = 'in')


ggplot(Thyroid_corr_TWAS_2, aes(x = Module, y = correlation.TWAS)) +
  geom_col(aes(fill = Module), colour = 'grey70') +
  geom_hline(yintercept = quantile(Thyroid_corr_TWAS_2$correlation.TWAS, my.quantile, na.rm = T)[2], color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = quantile(Thyroid_corr_TWAS_2$correlation.TWAS, my.quantile, na.rm = T)[1], color = 'darkblue', linetype = 'dashed') +
  # labs(title = "Correlation between Hyperthyroidism and Hypothyroidism") +
  ylab("Spearman's rho") + 
  theme_bw() +
  scale_fill_manual(values = color_df_2$Char) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6)) +
  facet_wrap(~Cell_type, nrow = 2, scales = 'free_x')

ggsave('./Grid_anal_073124/Correlation_TWAS_hyper_GD.png', width = 10, height = 3.5, dpi = 300, units = 'in')

ggplot(Thyroid_corr_delta_1, aes(x = Module, y = correlation.delta)) +
  geom_col(aes(fill = Module), colour = 'grey70') +
  geom_hline(yintercept = quantile(Thyroid_corr_delta_1$correlation.delta, my.quantile, na.rm = T)[2], color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = quantile(Thyroid_corr_delta_1$correlation.delta, my.quantile, na.rm = T)[1], color = 'darkblue', linetype = 'dashed') +
  # labs(title = "Correlation between Hyperthyroidism and Hypothyroidism") +
  ylab("Spearman's rho") + 
  theme_bw() +
  scale_fill_manual(values = color_df_1$Char) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6)) +
  facet_wrap(~Cell_type, nrow = 2, scales = 'free_x')

ggsave('./Grid_anal_073124/Correlation_delta_hypo_GD.png', width = 10, height = 3.5, dpi = 300, units = 'in')


ggplot(Thyroid_corr_delta_2, aes(x = Module, y = correlation.delta)) +
  geom_col(aes(fill = Module), colour = 'grey70') +
  geom_hline(yintercept = quantile(Thyroid_corr_delta_2$correlation.delta, my.quantile, na.rm = T)[2], color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = quantile(Thyroid_corr_delta_2$correlation.delta, my.quantile, na.rm = T)[1], color = 'darkblue', linetype = 'dashed') +
  # labs(title = "Correlation between Hyperthyroidism and Hypothyroidism") +
  ylab("Spearman's rho") + 
  theme_bw() +
  scale_fill_manual(values = color_df_2$Char) +
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6)) +
  facet_wrap(~Cell_type, nrow = 2, scales = 'free_x')

ggsave('./Grid_anal_073124/Correlation_delta_hyper_GD.png', width = 10, height = 3.5, dpi = 300, units = 'in')

# Thyroid_corr_delta <- Thyroid_corr_delta %>%
#   as.data.frame()

Thyroid_corr_delta_1 <- Thyroid_corr_delta_1 %>%
  dplyr::mutate(., module.interest = ifelse(correlation.delta > quantile(Thyroid_corr_delta_1$correlation.delta, my.quantile)[2],
                                            'Positive', ifelse(correlation.delta < quantile(Thyroid_corr_delta_1$correlation.delta, my.quantile)[1], 
                                                               'Negative', 
                                                               'None')))

Thyroid_corr_1 <- Thyroid_corr_1 %>%
  dplyr::mutate(., module.interest = ifelse(correlation.propagation > quantile(Thyroid_corr_1$correlation.propagation, my.quantile)[2],
                                            'Positive', ifelse(correlation.propagation < quantile(Thyroid_corr_1$correlation.propagation, my.quantile)[1], 
                                                               'Negative', 
                                                               'None')))


Thyroid_corr_delta_2 <- Thyroid_corr_delta_2 %>%
  dplyr::mutate(., module.interest = ifelse(correlation.delta > quantile(Thyroid_corr_delta_2$correlation.delta, my.quantile)[2],
                                            'Positive', ifelse(correlation.delta < quantile(Thyroid_corr_delta_2$correlation.delta, my.quantile)[1], 
                                                               'Negative', 
                                                               'None')))

Thyroid_corr_2 <- Thyroid_corr_2 %>%
  dplyr::mutate(., module.interest = ifelse(correlation.propagation > quantile(Thyroid_corr_2$correlation.propagation, my.quantile)[2],
                                            'Positive', ifelse(correlation.propagation < quantile(Thyroid_corr_2$correlation.propagation, my.quantile)[1], 
                                                               'Negative', 
                                                               'None')))

Thyroid_delta_modules_hypo_GD <- Thyroid_corr_delta_1 %>%
  dplyr::filter(.,module.interest != 'None')

Thyroid_delta_modules_hyper_GD <- Thyroid_corr_delta_2 %>%
  dplyr::filter(.,module.interest != 'None')



## Test plot of MOI based on delta value.

total_counts_1 <- Thyroid_corr_delta_1 %>%
  group_by(Cell_type) %>%
  summarise(total_count = n())

# Calculate the proportion of each module of interest within each Cell_type

proportions_1 <- Thyroid_corr_delta_1 %>%
  left_join(total_counts_1, by = "Cell_type") %>%
  group_by(Cell_type, module.interest) %>%
  summarise(proportion = n() / total_count) %>%
  mutate(label = paste0(round(proportion * 100, 1), "%")) %>% 
  unique

proportions_1$module.interest <- proportions_1$module.interest %>%
  factor(., levels = c('Positive', 'None', 'Negative'))

proportions_1 <- proportions_1 %>%
  mutate(., loc = if_else(module.interest == 'Positive',
                          0.005,
                          if_else(module.interest == 'Negative',
                                  0.972,
                                  0.55)
  ))


# Plot the bar plot

ggplot(Thyroid_corr_delta_1, aes(x = Cell_type, fill = module.interest)) +
  geom_bar(position = 'fill') +
  geom_text(data = proportions_1, aes(label = label, y = loc), 
            vjust = -0.5,
            size = 1.7) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('Relative Proportion') +
  labs(fill = 'Module Correlation') +
  scale_fill_manual(values = c('lightpink2', 'linen', 'lightsteelblue4')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))

ggsave(filename = './Grid_anal_073124/Corr_proportion_delta_barplot_GD_hypo.png', width = 5, height = 4, dpi = 300,
       units = 'in')



total_counts_2 <- Thyroid_corr_delta_2 %>%
  group_by(Cell_type) %>%
  summarise(total_count = n())

# Calculate the proportion of each module.interest within each Cell_type
proportions_2 <- Thyroid_corr_delta_2 %>%
  left_join(total_counts_2, by = "Cell_type") %>%
  group_by(Cell_type, module.interest) %>%
  summarise(proportion = n() / total_count) %>%
  mutate(label = paste0(round(proportion * 100, 1), "%")) %>% 
  unique

proportions_2$module.interest <- proportions_2$module.interest %>%
  factor(., levels = c('Positive', 'None', 'Negative'))

proportions_2 <- proportions_2 %>%
  mutate(., loc = if_else(module.interest == 'Positive',
                          0.012,
                          if_else(module.interest == 'Negative',
                                  0.95,
                                  0.38)
  ))


# Plot the bar plot

ggplot(Thyroid_corr_delta_2, aes(x = Cell_type, fill = module.interest)) +
  geom_bar(position = 'fill') +
  geom_text(data = proportions_2, aes(label = label, y = loc), 
            vjust = -0.5,
            size = 1.7) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('Relative Proportion') +
  labs(fill = 'Module Correlation') +
  scale_fill_manual(values = c('lightpink2', 'linen', 'lightsteelblue4')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))

ggsave(filename = './Grid_anal_073124/Corr_proportion_delta_barplot_GD_hyper.png', width = 5, height = 4, dpi = 300,
       units = 'in')



Thyroid_modules_hypo_GD <- Thyroid_corr_1 %>%
  dplyr::filter(.,module.interest != 'None')

Thyroid_modules_hyper_GD <- Thyroid_corr_2 %>%
  dplyr::filter(.,module.interest != 'None')



## Test plot of MOI based on Final.Heat value.

total_counts_1 <- Thyroid_corr_1 %>%
  group_by(Cell_type) %>%
  summarise(total_count = n())

# Calculate the proportion of each module of interest within each Cell_type

proportions_1 <- Thyroid_corr_1 %>%
  left_join(total_counts_1, by = "Cell_type") %>%
  group_by(Cell_type, module.interest) %>%
  summarise(proportion = n() / total_count) %>%
  mutate(label = paste0(round(proportion * 100, 1), "%")) %>% 
  unique

proportions_1$module.interest <- proportions_1$module.interest %>%
  factor(., levels = c('Positive', 'None', 'Negative'))

proportions_1 <- proportions_1 %>%
  mutate(., loc = if_else(module.interest == 'Positive',
                          0.004,
                          if_else(module.interest == 'Negative',
                                  0.976,
                                  0.55)
  ))


# Plot the bar plot

ggplot(Thyroid_corr_1, aes(x = Cell_type, fill = module.interest)) +
  geom_bar(position = 'fill') +
  geom_text(data = proportions_1, aes(label = label, y = loc), 
            vjust = -0.5,
            size = 1.7) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('Relative Proportion') +
  labs(fill = 'Module Correlation') +
  scale_fill_manual(values = c('lightpink2', 'linen', 'lightsteelblue4')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))

ggsave(filename = './Grid_anal_073124/Corr_proportion_FH_barplot_GD_hypo.png', width = 5, height = 4, dpi = 300,
       units = 'in')



total_counts_2 <- Thyroid_corr_2 %>%
  group_by(Cell_type) %>%
  summarise(total_count = n())

# Calculate the proportion of each module.interest within each Cell_type
proportions_2 <- Thyroid_corr_2 %>%
  left_join(total_counts_2, by = "Cell_type") %>%
  group_by(Cell_type, module.interest) %>%
  summarise(proportion = n() / total_count) %>%
  mutate(label = paste0(round(proportion * 100, 1), "%")) %>% 
  unique

proportions_2$module.interest <- proportions_2$module.interest %>%
  factor(., levels = c('Positive', 'None', 'Negative'))

proportions_2 <- proportions_2 %>%
  mutate(., loc = if_else(module.interest == 'Positive',
                          0.004,
                          if_else(module.interest == 'Negative',
                                  0.959,
                                  0.39)
  ))


# Plot the bar plot

ggplot(Thyroid_corr_2, aes(x = Cell_type, fill = module.interest)) +
  geom_bar(position = 'fill') +
  geom_text(data = proportions_2, aes(label = label, y = loc), 
            vjust = -0.5,
            size = 1.7) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('Relative Proportion') +
  labs(fill = 'Module Correlation') +
  scale_fill_manual(values = c('lightpink2', 'linen', 'lightsteelblue4')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))

ggsave(filename = './Grid_anal_073124/Corr_proportion_FH_barplot_GD_hyper.png', width = 5, height = 4, dpi = 300,
       units = 'in')


#### Combining module prioritizations with both differential value & Final heat

Thyroid_together_1 <- Thyroid_corr_1 %>%
  inner_join(., Thyroid_corr_delta_1, by = c('Module', 'Cell_type')) %>%
  dplyr::mutate(., module.interest = if_else(module.interest.x == 'Negative' & module.interest.y == 'Negative', 'Negative',
                                             if_else(module.interest.x == 'Positive' & module.interest.y == 'Positive', 'Positive',
                                                     'None')))

Thyroid_together_2 <- Thyroid_corr_2 %>%
  inner_join(., Thyroid_corr_delta_2, by = c('Module', 'Cell_type')) %>%
  dplyr::mutate(., module.interest = if_else(module.interest.x == 'Negative' & module.interest.y == 'Negative', 'Negative',
                                             if_else(module.interest.x == 'Positive' & module.interest.y == 'Positive', 'Positive',
                                                     'None')))



## Test plot of MOI based on Final.Heat value.

total_counts_1 <- Thyroid_together_1 %>%
  group_by(Cell_type) %>%
  summarise(total_count = n())

# Calculate the proportion of each module of interest within each Cell_type

proportions_1 <- Thyroid_together_1 %>%
  left_join(total_counts_1, by = "Cell_type") %>%
  group_by(Cell_type, module.interest) %>%
  summarise(proportion = n() / total_count) %>%
  mutate(label = paste0(round(proportion * 100, 1), "%")) %>% 
  unique

proportions_1$module.interest <- proportions_1$module.interest %>%
  factor(., levels = c('Positive', 'None', 'Negative'))

proportions_1 <- proportions_1 %>%
  mutate(., loc = if_else(module.interest == 'Positive',
                          -0.002,
                          if_else(module.interest == 'Negative',
                                  0.952,
                                  0.55)
  ))


# Plot the bar plot

ggplot(Thyroid_together_1, aes(x = Cell_type, fill = module.interest)) +
  geom_bar(position = 'fill') +
  geom_text(data = proportions_1, aes(label = label, y = loc), 
            vjust = -0.5,
            size = 1.7) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('Relative Proportion') +
  labs(fill = 'Module Correlation') +
  scale_fill_manual(values = c('lightpink2', 'linen', 'lightsteelblue4')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))

ggsave(filename = './Grid_anal_073124/Corr_proportion_together_barplot_GD_hypo.png', width = 5, height = 4, dpi = 300,
       units = 'in')


## For Hyper


total_counts_2 <- Thyroid_together_2 %>%
  group_by(Cell_type) %>%
  summarise(total_count = n())

# Calculate the proportion of each module.interest within each Cell_type
proportions_2 <- Thyroid_together_2 %>%
  left_join(total_counts_2, by = "Cell_type") %>%
  group_by(Cell_type, module.interest) %>%
  summarise(proportion = n() / total_count) %>%
  mutate(label = paste0(round(proportion * 100, 1), "%")) %>% 
  unique

proportions_2$module.interest <- proportions_2$module.interest %>%
  factor(., levels = c('Positive', 'None', 'Negative'))

proportions_2 <- proportions_2 %>%
  mutate(., loc = if_else(module.interest == 'Positive',
                          0.004,
                          if_else(module.interest == 'Negative',
                                  0.959,
                                  0.55)
  ))


# Plot the bar plot

ggplot(Thyroid_together_2, aes(x = Cell_type, fill = module.interest)) +
  geom_bar(position = 'fill') +
  geom_text(data = proportions_2, aes(label = label, y = loc), 
            vjust = -0.5,
            size = 1.7) +
  theme_bw() +
  xlab('Cell Type') +
  ylab('Relative Proportion') +
  labs(fill = 'Module Correlation') +
  scale_fill_manual(values = c('lightpink2', 'linen', 'lightsteelblue4')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))

ggsave(filename = './Grid_anal_073124/Corr_proportion_together_barplot_GD_hyper.png', width = 5, height = 4, dpi = 300,
       units = 'in')




## Narrowing prioritized genes in module of interests (MOI)
library(pheatmap)

library(ComplexHeatmap)

library(circlize)

head(MC_res)

# load('./073124_data_to_anal.Rdata')

load('./Sig_gene_per_method/Dat_with_significance.Rdata')

head(signif.tab)

signif.tab <- signif.tab %>%
  dplyr::filter(., method == 'mc') %>%
  dplyr::select(-Module_per_Cell)


## 1: hypo 2: hyper

hypo_modules <- Thyroid_together_1 %>%
  dplyr::filter(.,module.interest != 'None') %>%
  dplyr::mutate(., Phenotype = 'GD-Hypo', Cell_module = paste0(Cell_type, '_', Module))


hyper_modules <- Thyroid_together_2 %>%
  dplyr::filter(.,module.interest != 'None') %>%
  dplyr::mutate(., Phenotype = 'GD-Hyper', Cell_module = paste0(Cell_type, '_', Module))


share_modules <- intersect(hypo_modules$Cell_module, 
                           hyper_modules$Cell_module) 

Unmatched_modules <- union(setdiff(hypo_modules$Cell_module, hyper_modules$Cell_module),
                           setdiff(hyper_modules$Cell_module, hypo_modules$Cell_module))

All_module <- rbind(hypo_modules, hyper_modules) %>%
  dplyr::mutate(., is.shared = if_else(Cell_module %in% Unmatched_modules, 'Unmatced', 'Shared'))

Int_modules <- All_module %>%
  dplyr::filter(., is.shared == 'Shared')

share_modules <- Int_modules %>%
  dplyr::filter(., module.interest == 'Positive') %>%
  .$Cell_module %>% unique

Specific_modules <- Int_modules %>%
  dplyr::filter(., module.interest == 'Negative') %>%
  .$Cell_module %>% unique

signif.tab <- signif.tab %>%
  dplyr::mutate(., Cell_module = paste0(Cell_type, '_', Module))


Shared_mat <- signif.tab %>%
  dplyr::filter(., Cell_module %in% share_modules)


Specific_mat <- signif.tab %>%
  dplyr::filter(., Cell_module %in% Specific_modules)


Shared_mat$Cell_module <- factor(Shared_mat$Cell_module)

heat.shared.pre <- Shared_mat[order(Shared_mat$Cell_module), ]

shared.gene.info <- heat.shared.pre %>% dplyr::select(., gene, Cell_module, Cell_type, Module) %>% 
  dplyr::distinct()


hypo_shared <- MC_res %>%
  dplyr::filter(., Phenotype == 'hypo', gene %in% heat.shared.pre$gene) %>%
  semi_join(shared.gene.info, by = c('gene', 'Cell_type', 'Module'))

hyper_shared <- MC_res %>%
  dplyr::filter(., Phenotype == 'hyper', gene %in% heat.shared.pre$gene) %>%
  semi_join(shared.gene.info, by = c('gene', 'Cell_type', 'Module'))

heat.shared <- data.frame(hypo = hypo_shared$delta,
                          GD = heat.shared.pre$delta,
                          hyper = hyper_shared$delta,
                          gene = heat.shared.pre$gene,
                          Cell_module = heat.shared.pre$Cell_module)

heat.dir <- '/home/js/Thyroid_disorder/Grid_anal_073124/Shared_heat/'

# Ensure the directory exists
if (!dir.exists(heat.dir)) {
  dir.create(heat.dir)
}

# Close all open graphical devices at the start
while (dev.cur() > 1) dev.off()

# Iterate over unique cell modules
for (i in unique(as.character(heat.shared.pre$Cell_module))) {
  print(paste0("Processing Cell Module: ", i))  # Debugging message
  
  # Filter data for the current cell module
  tmp.frame <- heat.shared %>%
    dplyr::filter(Cell_module == i)
  
  rownames(tmp.frame) <- tmp.frame$gene
  tmp.frame <- tmp.frame[, -c(4, 5)]  # Exclude specific columns
  
  # Prepare annotation data
  tmp.anno <- shared.gene.info %>%
    dplyr::filter(Cell_module == i)
  
  write.table(tmp.anno$gene,
              file = paste0(heat.dir, '/', i, '.txt'),
              row.names = FALSE,
              col.names = FALSE,
              sep = '\n',
              quote = FALSE)
  
  # Create row annotation
  annotation_row <- data.frame(Module = as.character(tmp.anno$Module))
  rownames(annotation_row) <- rownames(tmp.frame)
  
  # Define annotation colors as a named list
  unique_modules <- unique(annotation_row$Module)
  annotation_colors <- list(Module = setNames(unique(annotation_row$Module), unique(annotation_row$Module)))
  
  
  # Define color gradient for the heatmap
  heatmap_colors <- colorRamp2(c(min(tmp.frame, na.rm = TRUE), 0, max(tmp.frame, na.rm = TRUE)), 
                               c("darkblue", "ivory", "firebrick3"))
  
  # Set output filename
  filename <- paste0(heat.dir, '/', i, '.png')
  print(paste0("Saving file: ", filename))  # Debugging message
  
  # Heatmap plot dimensions based on the number of rows
  height <- ifelse(nrow(tmp.frame) > 100, 0.01 * nrow(tmp.frame), 4.36)
  
  width <- ifelse(nrow(tmp.frame) > 100, 2.2, 2.5)
  
  border_col <- 'white'
  
  if(nrow(tmp.frame) > 100){
    border_lwd <- 0.001
  }else{
    border_lwd <- 0.8
  }
  
  
  # Open PNG device
  png(filename = filename, width = width, height = height, res = 300, units = 'in')
  
  # Create the heatmap
  heatmap <- Heatmap(as.matrix(tmp.frame),
                     name = "Delta",               # Legend title
                     col = heatmap_colors,              # Color palette
                     show_row_names = nrow(tmp.frame) <= 100,  # Show row names only for small matrices
                     show_column_names = TRUE,          # Always show column names
                     cluster_rows = FALSE,              # Do not cluster rows
                     cluster_columns = FALSE,           # Do not cluster columns
                     border = T,
                     rect_gp = gpar(col = border_col, lwd = border_lwd),
                     row_title = NA,               # Title for rows
                     column_title = NA,          # Title for columns
                     row_names_gp = gpar(fontsize = 4.5), # Font size for row names
                     column_names_gp = gpar(fontsize = 4.5), # Font size for column names
                     left_annotation = rowAnnotation(Module = annotation_row$Module,
                                                     col = annotation_colors,
                                                     border = T,
                                                     gp = gpar(col = border_col, lwd = border_lwd),
                                                     annotation_name_gp = gpar(fontsize = 4.5),
                                                     annotation_legend_param = list(labels_gp = gpar(fontsize = 4.5),
                                                                                    title_gp = gpar(fontsize = 4.5),
                                                                                    legend_height = unit(0.1, 'in'),
                                                                                    grid_height = unit(0.3, 'cm'),
                                                                                    grid_width = unit(0.3, 'cm'))
                                                     ),
                     heatmap_legend_param = list(direction = 'horizontal',
                                                 title_gp = gpar(fontsize = 4.5),
                                                 labels_gp = gpar(fontsize = 4.5),
                                                 legend_height = unit(0.1, 'in'),
                                                 legend_width = unit(1, 'in'),
                                                 grid_height = unit(0.2, 'cm'),
                                                 grid_width = unit(0.5, 'in')))  
  
  # Draw the heatmap
  draw(heatmap, heatmap_legend_side = "top")  # Legend below the heatmap
  
  # Close the PNG device
  dev.off()
}


# Next, visualizing specific modules.

Specific_mat$Cell_module <- factor(Specific_mat$Cell_module)

heat.specific.pre <- Specific_mat[order(Specific_mat$Cell_module), ]

shared.gene.info <- heat.specific.pre %>% dplyr::select(., gene, Cell_type, Cell_module, Module) %>% 
  dplyr::distinct()


hypo_specific <- MC_res %>%
  dplyr::filter(., Phenotype == 'hypo', gene %in% heat.specific.pre$gene) %>%
  semi_join(shared.gene.info, by = c('gene', 'Cell_type', 'Module'))

hyper_specific <- MC_res %>%
  dplyr::filter(., Phenotype == 'hyper', gene %in% heat.specific.pre$gene) %>%
  semi_join(shared.gene.info, by = c('gene', 'Cell_type', 'Module'))

heat.specific <- data.frame(hypo = hypo_specific$delta,
                            GD = heat.specific.pre$delta,
                            hyper = hyper_specific$delta,
                            gene = heat.specific.pre$gene,
                            Cell_module = heat.specific.pre$Cell_module)


heat.dir <- '/home/js/Thyroid_disorder/Grid_anal_073124/Specific_heat/'


if (!dir.exists(heat.dir)) {
  dir.create(heat.dir)
}

# Close all open graphical devices at the start
while (dev.cur() > 1) dev.off()

# Iterate over unique cell modules
for (i in unique(as.character(heat.specific.pre$Cell_module))) {
  print(paste0("Processing Cell Module: ", i))  # Debugging message
  
  # Filter data for the current cell module
  tmp.frame <- heat.specific %>%
    dplyr::filter(Cell_module == i)
  
  rownames(tmp.frame) <- tmp.frame$gene
  tmp.frame <- tmp.frame[, -c(4, 5)]  # Exclude specific columns
  
  # Prepare annotation data
  tmp.anno <- shared.gene.info %>%
    dplyr::filter(Cell_module == i)
  
  write.table(tmp.anno$gene,
              file = paste0(heat.dir, '/', i, '.txt'),
              row.names = FALSE,
              col.names = FALSE,
              sep = '\n',
              quote = FALSE)
  
  # Create row annotation
  annotation_row <- data.frame(Module = as.character(tmp.anno$Module))
  rownames(annotation_row) <- rownames(tmp.frame)
  
  # Define annotation colors as a named list
  unique_modules <- unique(annotation_row$Module)
  annotation_colors <- list(Module = setNames(unique(annotation_row$Module), unique(annotation_row$Module)))
  
  
  # Define color gradient for the heatmap
  heatmap_colors <- colorRamp2(c(min(tmp.frame, na.rm = TRUE), 0, max(tmp.frame, na.rm = TRUE)), 
                               c("darkblue", "ivory", "firebrick3"))
  
  # Set output filename
  filename <- paste0(heat.dir, '/', i, '.png')
  print(paste0("Saving file: ", filename))  # Debugging message
  
  # Heatmap plot dimensions based on the number of rows
  height <- ifelse(nrow(tmp.frame) > 100, 0.01 * nrow(tmp.frame), 4.36)
  
  # Open PNG device
  png(filename = filename, width = 2.5, height = height, res = 300, units = 'in')
  
  # Create the heatmap
  heatmap <- Heatmap(as.matrix(tmp.frame),
                     name = "Delta",               # Legend title
                     col = heatmap_colors,              # Color palette
                     show_row_names = nrow(tmp.frame) <= 100,  # Show row names only for small matrices
                     show_column_names = TRUE,          # Always show column names
                     cluster_rows = FALSE,              # Do not cluster rows
                     cluster_columns = FALSE,           # Do not cluster columns
                     border = T,
                     rect_gp = gpar(col = 'white', lwd = 0.8),
                     row_title = NA,               # Title for rows
                     column_title = NA,          # Title for columns
                     row_names_gp = gpar(fontsize = 3.5), # Font size for row names
                     column_names_gp = gpar(fontsize = 4.5), # Font size for column names
                     left_annotation = rowAnnotation(Module = annotation_row$Module,
                                                     col = annotation_colors,
                                                     border = T,
                                                     gp = gpar(col = border_col, lwd = border_lwd),
                                                     annotation_name_gp = gpar(fontsize = 4.5),
                                                     annotation_legend_param = list(labels_gp = gpar(fontsize = 4.5),
                                                                                    title_gp = gpar(fontsize = 4.5),
                                                                                    legend_height = unit(0.1, 'in'),
                                                                                    grid_height = unit(0.3, 'cm'),
                                                                                    grid_width = unit(0.3, 'cm'))
                     ),
                     heatmap_legend_param = list(direction = 'horizontal',
                                                 title_gp = gpar(fontsize = 4.5),
                                                 labels_gp = gpar(fontsize = 4.5),
                                                 legend_height = unit(0.1, 'in'),
                                                 legend_width = unit(1, 'in'),
                                                 grid_height = unit(0.2, 'cm'),
                                                 grid_width = unit(0.5, 'in')))  
  
  # Draw the heatmap
  draw(heatmap, heatmap_legend_side = "top")  # Legend below the heatmap
  
  # Close the PNG device
  dev.off()
}


save(file = './Grid_anal_073124/Values.RData', heat.shared.pre, heat.specific.pre)

# Trying to visualize the results with dot plots.
my.quantile


GD_hypo <- merge(Thyroid_corr_1, Thyroid_corr_delta_1, by = c('Module', 'Cell_type'))

GD_hypo <- GD_hypo %>%
  dplyr::mutate(., Cell_module = paste0(Cell_type, ' (', Module, ')')) %>%
  dplyr::mutate(., Comb.signif = if_else(module.interest.x == 'Positive' & module.interest.y == 'Positive', 'Positive',
                                         if_else(module.interest.x == 'Negative' & module.interest.y == 'Negative', 'Negative', 
                                                 'None')))

ggplot() +
  geom_point(data = dplyr::filter(GD_hypo, Comb.signif == 'Positive'), 
             aes(x = correlation.propagation, y = correlation.delta), colour = 'purple3', size = 3, alpha = 0.9) +
  geom_point(data = dplyr::filter(GD_hypo, Comb.signif == 'Negative'), 
             aes(x = correlation.propagation, y = correlation.delta), colour = 'orange2', size = 3, alpha = 0.9) +
  geom_point(data = dplyr::filter(GD_hypo, Comb.signif == 'None'), 
             aes(x = correlation.propagation, y = correlation.delta), colour = 'grey60', alpha = 0.6) +
  geom_hline(yintercept = quantile(GD_hypo$correlation.delta, my.quantile, na.rm = T)[2], 
             color = 'red', linetype = 'dashed', alpha = 0.6) +
  geom_hline(yintercept = quantile(GD_hypo$correlation.delta, my.quantile, na.rm = T)[1], 
             color = 'darkblue', linetype = 'dashed', alpha = 0.6) +
  geom_vline(xintercept = quantile(GD_hypo$correlation.propagation, my.quantile, na.rm = T)[2], 
             color = 'red', linetype = 'dashed', alpha = 0.6) +
  geom_vline(xintercept = quantile(GD_hypo$correlation.propagation, my.quantile, na.rm = T)[1], 
             color = 'darkblue', linetype = 'dashed', alpha = 0.6) +
  theme_bw() +
  ylab('Delta') +
  xlab('Final Heat') + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6))

ggsave(filename = './Figs_for_submission/GD_hypo_module_corr.png', width = 3.5, height = 3.3, units = 'in', dpi = 300)
  

GD_hyper <- merge(Thyroid_corr_2, Thyroid_corr_delta_2, by = c('Module', 'Cell_type'))

GD_hyper <- GD_hyper %>%
  dplyr::mutate(., Cell_module = paste0(Cell_type, ' (', Module, ')')) %>%
  dplyr::mutate(., Comb.signif = if_else(module.interest.x == 'Positive' & module.interest.y == 'Positive', 'Positive',
                                         if_else(module.interest.x == 'Negative' & module.interest.y == 'Negative', 'Negative', 
                                                 'None')))

ggplot() +
  geom_point(data = dplyr::filter(GD_hyper, Comb.signif == 'Positive'), 
             aes(x = correlation.propagation, y = correlation.delta), colour = 'purple3', size = 3, alpha = 0.9) +
  geom_point(data = dplyr::filter(GD_hyper, Comb.signif == 'Negative'), 
             aes(x = correlation.propagation, y = correlation.delta), colour = 'orange2', size = 3, alpha = 0.9) +
  geom_point(data = dplyr::filter(GD_hyper, Comb.signif == 'None'), 
             aes(x = correlation.propagation, y = correlation.delta), colour = 'grey60', alpha = 0.6) +
  geom_hline(yintercept = quantile(GD_hyper$correlation.delta, my.quantile, na.rm = T)[2], 
             color = 'red', linetype = 'dashed', alpha = 0.6) +
  geom_hline(yintercept = quantile(GD_hyper$correlation.delta, my.quantile, na.rm = T)[1], 
             color = 'darkblue', linetype = 'dashed', alpha = 0.6) +
  geom_vline(xintercept = quantile(GD_hyper$correlation.propagation, my.quantile, na.rm = T)[2], 
             color = 'red', linetype = 'dashed', alpha = 0.6) +
  geom_vline(xintercept = quantile(GD_hyper$correlation.propagation, my.quantile, na.rm = T)[1], 
             color = 'darkblue', linetype = 'dashed', alpha = 0.6) +
  theme_bw() +
  ylab('Delta') +
  xlab('Final Heat') + 
  theme(legend.position = 'none', axis.text.y = element_text(size = 7),    
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 6))

ggsave(filename = './Figs_for_submission/GD_hyper_module_corr.png', width = 3.5, height = 3.3, units = 'in', dpi = 300)

## Now, let's make an UpSet plot with them!!        
library(ComplexHeatmap)

GD_hypo_pos <- GD_hypo %>%
  dplyr::filter(., Comb.signif == 'Positive') %>%
  .$Cell_module
  
GD_hypo_neg <- GD_hypo %>%
  dplyr::filter(., Comb.signif == 'Negative') %>%
  .$Cell_module

GD_hypo_non <- GD_hypo %>%
  dplyr::filter(., Comb.signif == 'None') %>%
  .$Cell_module


GD_hyper_pos <- GD_hyper %>%
  dplyr::filter(., Comb.signif == 'Positive') %>%
  .$Cell_module

GD_hyper_neg <- GD_hyper %>%
  dplyr::filter(., Comb.signif == 'Negative') %>%
  .$Cell_module

GD_hyper_non <- GD_hyper %>%
  dplyr::filter(., Comb.signif == 'None') %>%
  .$Cell_module


modules_lt <- list(`GD-hypothyroidism \n(Positive)` = GD_hypo_pos,
                   `GD-hypothyroidism \n(Negative)` = GD_hypo_neg,
                   `GD-hypothyroidism \n(None)` = GD_hypo_non,
                   `GD-hyperthyroidism \n(Postivie)` = GD_hyper_pos,
                   `GD-hyperthyroidism \n(Negative)` = GD_hyper_neg,
                   `GD-hypethyroidism \n(None)` = GD_hyper_non) %>%
  make_comb_mat()

png(filename = './Figs_for_submission/UpSet_Modules_0119.png', width = 9.5, height = 4.0, res = 300, units = 'in')

UpSet(modules_lt, lwd = 4, set_order = rownames(modules_lt),
      top_annotation = upset_top_annotation(modules_lt, add_numbers = TRUE),
      right_annotation = upset_right_annotation(modules_lt, add_numbers = TRUE, annotation_name_side = "top",
                                                axis_param = list(side = "top")),
      comb_col = c('purple3', 'darkseagreen4', 'orange2', 'darkseagreen4', 'darkseagreen4', 'darkseagreen4', 
                   'darkseagreen4', 'darkseagreen4'))

dev.off()


