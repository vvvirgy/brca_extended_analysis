library(tidyverse)
source('R/utils.R')
muts_selected = readRDS('data/r_obj/metabric/data_muts/muts_selected.rds')
brca_drivers = readRDS('../analysis/r_obj/brca_drivers.rds')

muts_selected_drivers = muts_selected %>% 
  filter(Hugo_Symbol %in% brca_drivers) %>% 
  select(Hugo_Symbol, HGVSp, Chromosome, Start_Position, End_Position, HGVSp_Short) %>% 
  distinct()

saveRDS(muts_selected_drivers, 'data/muts_selected_drivers.rds')
write.table(muts_selected_drivers, 'data/muts_selected_drivers.csv', sep = ',', quote = F, col.names = T, row.names = F)

metadata_selected = readRDS('data/r_obj/metabric/metadata/metadata_selected_samples.rds')
metadata_patient = readRDS('data/r_obj/metabric/metadata/metadata_selected_patients.rds')

colnames(metadata_selected)
metadata_patient = metadata_patient %>% 
  mutate(INFERRED_MENOPAUSAL_STATE = ifelse(INFERRED_MENOPAUSAL_STATE == '', 'Unknown', INFERRED_MENOPAUSAL_STATE))

all_metadata = full_join(metadata_patient, metadata_selected, by = join_by('PATIENT_ID' == 'patient_id'))

muts_selected = muts_selected %>% 
  mutate(is_driver = ifelse(Hugo_Symbol %in% brca_drivers, TRUE, FALSE)) %>% 
  mutate(is_driver = ifelse(Variant_Classification %in% c('Silent', 'Intron'), FALSE, is_driver))

# add metadata of menopause

menopause_status = all_metadata %>% 
  select(sample_id, INFERRED_MENOPAUSAL_STATE)

muts_selected = muts_selected %>% 
  full_join(., menopause_status, by = join_by('Tumor_Sample_Barcode' == 'sample_id')) %>% 
  mutate(INFERRED_MENOPAUSAL_STATE = factor(INFERRED_MENOPAUSAL_STATE, levels = c('Pre', 'Post', 'Unknown'))) 
  
# get drivers
drivers = muts_selected %>% 
  filter(is_driver)

# plot drivers distribution by menopause status
drivers_by_meno_prop = get_props(drivers, group = 'INFERRED_MENOPAUSAL_STATE') 
drivers_by_menopause = plot_prop_samples(drivers_by_meno_prop, colors = menopause_cols)
saveRDS(drivers_by_menopause, 'results/metabric/mutations/drivers_by_menopause_plot.rds')
ggsave('results/metabric/mutations/drivers_percentage_groups.pdf', plot = drivers_by_menopause, width = 6, height = 8)

# drivers_by_menopause = drivers %>% 
#   # mutate(INFERRED_MENOPAUSAL_STATE = factor(INFERRED_MENOPAUSAL_STATE, levels = c('Pre', 'Post', 'Unknown'))) %>% 
#   group_by(INFERRED_MENOPAUSAL_STATE) %>% 
#   mutate(tot_samples = n()) %>%
#   group_by(INFERRED_MENOPAUSAL_STATE, Hugo_Symbol) %>% 
#   mutate(n = n(), prop = n/tot_samples) %>% 
#   ggplot(aes(y = Hugo_Symbol, fill = INFERRED_MENOPAUSAL_STATE, x = prop)) + 
#   geom_bar(stat = 'identity', position = 'dodge') + 
#   theme_bw() + 
#   scale_x_continuous(labels = scales::percent) + 
#   scale_fill_manual(values = menopause_cols) + 
#   labs(x = 'Percentage of samples', y = 'Gene symbols') + 
#   theme(legend.position = 'right') + 
#   guides(fill = guide_legend(title = 'Menopause status', nrow = 3))

drivers_by_menopause_2_df = get_props(drivers, group = 'Hugo_Symbol') 
drivers_by_menopause_2_plot = plot_prop_genes(drivers_by_menopause_2_df, colors = menopause_cols)
saveRDS(drivers_by_menopause_2_plot, 'results/metabric/mutations/drivers_by_menopause_2_plot.rds')
ggsave(plot = drivers_by_menopause_2_plot, 'results/metabric/mutations/drivers_percentage_groups_2.pdf', width = 6, height = 8)

# look at all genes without selecting drivers

genes_by_meno_props = get_props(muts_selected, group = 'INFERRED_MENOPAUSAL_STATE') %>% 
  filter(!is.na(Hugo_Symbol))
genes_plot = genes_by_meno_props %>% 
  group_by(INFERRED_MENOPAUSAL_STATE) %>% 
  slice_max(prop, n = 50) %>% 
  pull(Hugo_Symbol) %>% unique
genes_by_meno_props_reduced = genes_by_meno_props %>% 
  filter(Hugo_Symbol %in% genes_plot)
plot_prop_samples(genes_by_meno_props_reduced, colors = menopause_cols)

muts_selected %>% 
  get_props(., 'Hugo_Symbol') %>% 
  plot_prop_genes(., colors = menopause_cols)


# drivers_by_menopause_2_plot = drivers_by_menopause_2_df %>% 
#   ggplot(aes(y = Hugo_Symbol, fill = INFERRED_MENOPAUSAL_STATE, x = prop)) + 
#   geom_bar(stat = 'identity') + 
#   theme_bw() + 
#   scale_x_continuous(labels = scales::percent) +
#   scale_fill_manual(values = menopause_cols) + 
#   labs(x = 'Percentage of samples mutated', y = 'Gene symbols') + 
#   theme(legend.position = 'right') + 
#   guides(fill = guide_legend(title = 'Menopause status', nrow = 3))




# ggsave('results/metabric/mutations/drivers_percentage_groups.pdf', plot = drivers_by_menopause, width = 6, height = 8)


