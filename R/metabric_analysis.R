rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')
setwd('/orfeo/cephfs/scratch/area/vgazziero/CDSlab/brca_cro_prj/brca_extended_analysis')
library(tidyverse)
library(ggsurvfit)
library(survival)
library(wesanderson)
library(survminer)
source('R/utils.R')

# Preparation ---------

muts_selected = readRDS('data/r_obj/metabric/data_muts/muts_selected.rds')
brca_drivers = readRDS('../analysis/r_obj/brca_drivers.rds')

muts_selected_drivers = muts_selected %>% 
  filter(Hugo_Symbol %in% brca_drivers) %>% 
  dplyr::select(Hugo_Symbol, HGVSp, Chromosome, Start_Position, End_Position, starts_with('HGVS'), Consequence, Variant_Classification) %>% 
  distinct() 
# saveRDS(muts_selected_drivers, 'data/muts_selected_drivers.rds')
# write.table(muts_selected_drivers, 'data/muts_selected_drivers.csv', sep = ',', quote = F, col.names = T, row.names = F)

gata3_tp53 = muts_selected_drivers %>% 
  filter(Hugo_Symbol %in% c('TP53', 'GATA3')) %>% 
  distinct()

# saveRDS(gata3_tp53, 'data/gata3_tp53.rds')
# write.table(gata3_tp53, 'data/gata3_tp53.csv', sep = ',', quote = F, col.names = T, row.names = F)

metadata_selected = readRDS('data/r_obj/metabric/metadata/metadata_selected_samples.rds')
metadata_patient = readRDS('data/r_obj/metabric/metadata/metadata_selected_patients.rds')

# classify the menopausal status based on age
metadata_patient = metadata_patient %>% 
  mutate(Inferred_menopause_status = case_when(AGE_AT_DIAGNOSIS >= 55 ~ 'Post', 
                                               AGE_AT_DIAGNOSIS <= 55 ~ 'Pre', 
                                               .default = 'Unknown'))

all_metadata = full_join(metadata_patient, metadata_selected, by = join_by('PATIENT_ID' == 'patient_id'))
all_metadata = all_metadata %>% 
  mutate(Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post')))

saveRDS(all_metadata, 'data/r_obj/metabric/metadata/all_metadata.rds')

muts_selected = muts_selected %>% 
  mutate(is_driver = ifelse(Hugo_Symbol %in% brca_drivers, TRUE, FALSE)) %>% 
  mutate(is_driver = ifelse(Variant_Classification %in% c('Silent', 'Intron'), FALSE, is_driver))

# add metadata of menopause
menopause_status = all_metadata %>% 
  dplyr::select(sample_id, Inferred_menopause_status)

muts_selected = muts_selected %>% 
  full_join(., menopause_status, by = join_by('Tumor_Sample_Barcode' == 'sample_id')) %>% 
  mutate(Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post', 'Unknown'))) 
  
# SNV plots
## plots of drivers ----------  

# get drivers
drivers = muts_selected %>% 
  filter(is_driver)

# plot drivers distribution by menopause status
drivers_by_meno_prop = get_props(drivers, group = 'Inferred_menopause_status') 
drivers_by_menopause = plot_prop_samples(drivers_by_meno_prop, colors = menopause_cols)
saveRDS(drivers_by_menopause, 'results/metabric/mutations/drivers_by_menopause_plot.rds')
ggsave('results/metabric/mutations/drivers_percentage_groups.pdf', plot = drivers_by_menopause, width = 6, height = 8)

drivers_by_menopause_2_df = get_props(drivers, group = 'Hugo_Symbol') 
drivers_by_menopause_2_plot = plot_prop_genes(drivers_by_menopause_2_df, colors = menopause_cols)
saveRDS(drivers_by_menopause_2_plot, 'results/metabric/mutations/drivers_by_menopause_2_plot.rds')
ggsave(plot = drivers_by_menopause_2_plot, 'results/metabric/mutations/drivers_percentage_groups_2.pdf', width = 6, height = 8)

## look at all genes without selecting drivers ----
genes_by_meno_props = get_props(muts_selected, group = 'Inferred_menopause_status') %>% 
  filter(!is.na(Hugo_Symbol))
genes_plot = genes_by_meno_props %>% 
  group_by(Inferred_menopause_status) %>% 
  slice_max(prop, n = 50) %>% 
  pull(Hugo_Symbol) %>% unique
genes_by_meno_props_reduced = genes_by_meno_props %>% 
  filter(Hugo_Symbol %in% genes_plot)
plot_prop_samples(genes_by_meno_props_reduced, colors = menopause_cols)

muts_selected %>% 
  get_props(., 'Hugo_Symbol') %>% 
  plot_prop_genes(., colors = menopause_cols)


## taking a look to tmb distribution ----
tmb_boxplot = all_metadata %>% 
  rename(TMB = `TMB (nonsynonymous)`) %>% 
  mutate(TMB = as.numeric(TMB)) %>% 
  ggplot(aes(y = TMB, x = Inferred_menopause_status, 
             # color = INFERRED_MENOPAUSAL_STATE, 
             fill = Inferred_menopause_status)) +
  # geom_jitter(alpha = 0.1, inherit.aes = T, show.legend = F) +
  # geom_boxplot(outliers = F, alpha = 0.8) + 
  geom_violin(alpha = .8) + 
  ggpubr::stat_compare_means(comparisons = list(c('Pre', "Post"))) + 
  theme_bw() + 
  scale_fill_manual(values = menopause_cols) + 
  labs(x = 'Inferred menopausal state') + 
  theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(title = 'Inferred menopausal state'))
ggsave('results/metabric/mutations/TMB_boxplot_comparison.png', plot = tmb_boxplot, width = 6, height = 6)
saveRDS(tmb_boxplot, 'results/metabric/tmb_boxplot.rds')  


# SURVIVAL ANALYSIS --------------- 
# trying to perform some survival analysis

tp53_mutated = muts_selected %>% 
  filter(Hugo_Symbol == 'TP53') %>%
  group_by(Tumor_Sample_Barcode) %>% 
  mutate(multihit = ifelse(n() > 1, TRUE, FALSE))

tp53_mutated_samples = tp53_mutated$Tumor_Sample_Barcode %>% unique 

gata3_mutated = muts_selected %>% 
  filter(Hugo_Symbol == 'GATA3') %>%
  group_by(Tumor_Sample_Barcode) %>% 
  mutate(multihit = ifelse(n() > 1, TRUE, FALSE))

gata3_mutated_samples = gata3_mutated$Tumor_Sample_Barcode %>% unique 

# including also co-mutation of pik3ca and tp53
pik3ca_mutated = muts_selected %>% 
  filter(Hugo_Symbol == 'PIK3CA') %>%
  group_by(Tumor_Sample_Barcode) %>% 
  mutate(multihit = ifelse(n() > 1, TRUE, FALSE))

pik3ca_mutated_samples = pik3ca_mutated$Tumor_Sample_Barcode %>% unique


all_metadata = all_metadata %>% 
  mutate(TP53 = ifelse(PATIENT_ID %in% tp53_mutated_samples, 'Mutated', 'Wild-type'), 
         GATA3 = ifelse(PATIENT_ID %in% gata3_mutated_samples, 'Mutated', 'Wild-type'), 
         PIK3CA = ifelse(PATIENT_ID %in% pik3ca_mutated_samples, 'Mutated', 'Wild-type')
         ) %>% 
  mutate(comut_TP53_PIK3CA = ifelse(TP53 == 'Mutated' & PIK3CA == 'Mutated', 'co-mutated', 'absence of co-mutation'))

## by gene mutation status ----

survival_metadata = all_metadata %>% 
  dplyr::select(Inferred_menopause_status, OS_MONTHS, OS_STATUS, VITAL_STATUS, TP53, GATA3, PIK3CA, comut_TP53_PIK3CA) %>% 
  tidyr::separate(OS_STATUS, into = c('os_status', 'status'), sep = ':', convert = T) %>% 
  filter(VITAL_STATUS != "Died of Other Causes") %>% 
  filter(VITAL_STATUS != "") %>% 
  filter(os_status != '') %>% 
  mutate(os_years = OS_MONTHS/12) %>% 
  mutate(os_years = as.numeric(os_years))

survival_metadata = survival_metadata %>% 
  mutate(
    TP53 = factor(TP53, levels = c('Wild-type', 'Mutated')), 
    GATA3 = factor(GATA3, levels = c('Wild-type', 'Mutated')), 
    PIK3CA = factor(PIK3CA, levels = c('Wild-type', 'Mutated')), 
    comut_TP53_PIK3CA = factor(comut_TP53_PIK3CA, levels = c('absence of co-mutation', 'co-mutated')), 
    Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post'))
  )

os_analysis = lapply(c('TP53', 'GATA3', 'comut_TP53_PIK3CA'), function(g) {
  survival_analysis(survival_metadata, 
                    time = 'os_years', 
                    event = 'os_status', 
                    what = paste(g, '+ Inferred_menopause_status'), 
                    gene = g)
})

names(os_analysis) = c('TP53', 'GATA3', 'comut_TP53_PIK3CA')
saveRDS(os_analysis, 'results/metabric/mutations/gata3_tp53_survival.rds')

os_analysis$comut_TP53_PIK3CA$hr_plot
os_analysis$comut_TP53_PIK3CA$coxph %>% summary
  
os_analysis$TP53$hr_plot
os_analysis$GATA3$hr_plot


## by gene mutation status (indipendent pre and post) ----

survival_metadata = all_metadata %>% 
  dplyr::select(Inferred_menopause_status, OS_MONTHS, OS_STATUS, VITAL_STATUS, TP53, GATA3, PIK3CA, comut_TP53_PIK3CA) %>% 
  tidyr::separate(OS_STATUS, into = c('os_status', 'status'), sep = ':', convert = T) %>% 
  filter(VITAL_STATUS != "Died of Other Causes") %>% 
  filter(VITAL_STATUS != "") %>% 
  filter(os_status != '') %>% 
  mutate(os_years = OS_MONTHS/12) %>% 
  mutate(os_years = as.numeric(os_years)) %>% 
  mutate(os_status = case_when(
    
    (os_years >= 10 & os_status == 1) ~ 0, # censor the patients that were still alive 10y later
    .default = os_status
    
  ), 
  os_years = ifelse(os_years > 10, 10, os_years)) %>% 
  mutate(
    TP53 = factor(TP53, levels = c('Wild-type', 'Mutated')), 
    GATA3 = factor(GATA3, levels = c('Wild-type', 'Mutated')), 
    PIK3CA = factor(PIK3CA, levels = c('Wild-type', 'Mutated')), 
    comut_TP53_PIK3CA = factor(comut_TP53_PIK3CA, levels = c('absence of co-mutation', 'co-mutated')), 
    Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post'))
  ) %>% 
  # group_by(Inferred_menopause_status) %>% 
  base::split(.$Inferred_menopause_status) 

os_analysis = lapply(c('TP53', 'GATA3', 'comut_TP53_PIK3CA'), function(g) {
  lapply(survival_metadata, function(stat) {
    survival_analysis(stat, 
                      time = 'os_years', 
                      event = 'os_status', 
                      what = g, 
                      gene = g)
  })
})

names(os_analysis) = c('TP53', 'GATA3', 'comut_TP53_PIK3CA')
os_analysis$TP53$Pre$hr_plot
os_analysis$TP53$Post$hr_plot

os_analysis$GATA3$Pre$km_plot
os_analysis$GATA3$Post$km_plot
os_analysis$TP53$Post$km_plot
os_analysis$TP53$Pre$km_plot

saveRDS(os_analysis, 'results/metabric/mutations/gata3_tp53_survival_by_group.rds')
x = readRDS('results/metabric/mutations/gata3_tp53_survival_by_group.rds')
x$GATA3$Pre$plot
## mutation type of gata3 and tp53 ----

gata3_tp53_muts_classification = read.csv('data/gata3_tp53_domain_mutations_TP53.xlsx - mut_list_with_class.csv')

### GATA3 ----
gata3_cls = gata3_tp53_muts_classification %>% 
  filter(Hugo_Symbol == 'GATA3') %>% 
  mutate(Chromosome = as.character(Chromosome)) %>% 
  dplyr::select(-c(Consequence, Variant_Classification, X.1, X))

gata3_mutated_samples_cls = muts_selected %>% 
  full_join(., gata3_cls) %>% 
  filter(Hugo_Symbol == 'GATA3') %>%
  dplyr::select(Variant_Classification, Domain_Mutations, Tumor_Sample_Barcode) %>%
  rename(GATA3_mut_type = Domain_Mutations)

gata3_mutated_samples_cls = gata3_mutated_samples_cls %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  mutate(Variant_Classification = ifelse(n() > 1, 'multihit', Variant_Classification)) %>% 
  distinct() %>% 
  filter(n() == 1)


gata3_metadata = all_metadata %>% 
  dplyr::select(PATIENT_ID, Inferred_menopause_status, OS_MONTHS, OS_STATUS, VITAL_STATUS, TP53, GATA3, PIK3CA, comut_TP53_PIK3CA) %>% 
  tidyr::separate(OS_STATUS, into = c('os_status', 'status'), sep = ':', convert = T) %>% 
  filter(VITAL_STATUS != "Died of Other Causes") %>% 
  filter(VITAL_STATUS != "") %>% 
  filter(os_status != '') %>% 
  mutate(os_years = OS_MONTHS/12) %>% 
  mutate(os_years = as.numeric(os_years)) %>% 
  mutate(os_status = case_when(
    
    (os_years >= 10 & os_status == 1) ~ 0, # censor the patients that were still alive 10y later
    .default = os_status
    
  ), 
      os_years = ifelse(os_years > 10, 10, os_years)
  ) %>% 
  full_join(., gata3_mutated_samples_cls, by = join_by('PATIENT_ID' == 'Tumor_Sample_Barcode')) %>% 
  filter(!is.na(PATIENT_ID)) %>%
  filter(!is.na(os_status)) %>%
  filter(!is.na(VITAL_STATUS))

gata3_metadata = gata3_metadata %>% 
  mutate(GATA3_mut_type = case_when(
    (is.na(GATA3_mut_type) & GATA3 == 'Wild-type') ~ 'Wild-type/Silent mutation', 
    (is.na(GATA3_mut_type) & GATA3 == 'Mutated' & Variant_Classification == 'Silent') ~ 'Wild-type/Silent mutation', 
    .default = GATA3_mut_type
  ))

gata3_metadata = gata3_metadata %>% 
  mutate(GATA3_mut_type = factor(GATA3_mut_type, levels = c("Wild-type/Silent mutation", 'Other', 'ZnFn2_Retained', 'ZnFn2_Loss', 'ZnFn2_Modified'))) %>% 
  split(.$Inferred_menopause_status)

os_analysis_gata3_mut_type = lapply(c('GATA3_mut_type'), function(g) {
  lapply(gata3_metadata, function(stat) {
    survival_analysis(stat, 
                      time = 'os_years', 
                      event = 'os_status', 
                      what = g, 
                      gene = g)
  })
})

pre = os_analysis_gata3_mut_type[[1]]$Pre$km_plot + 
  ggtitle('Pre Menopause') + 
  labs(y = 'OS')
post = os_analysis_gata3_mut_type[[1]]$Post$km_plot + 
  ggtitle('Post Menopause') + 
  labs(y = 'OS')

pp = pre + post  
ggsave('results/metabric/mutations/gata3_mut_type_pre_post.png', width = 10, height = 6, plot = pp)
ggsave('results/metabric/mutations/gata3_mut_type_pre_post.pdf', width = 10, height = 6, plot = pp)

os_analysis_gata3_mut_type[[1]]$Post$hr_plot
os_analysis_gata3_mut_type[[1]]$Pre$hr_plot

# pdf('results/metabric/mutations/gata3_mut_type_pre_post.pdf', width = 10, height = 6)

pre = os_analysis_gata3_mut_type[[1]]$Pre$hr_plot + 
  ggtitle('Pre Menopause') + 
  labs(y = 'OS')
post = os_analysis_gata3_mut_type[[1]]$Post$hr_plot + 
  ggtitle('Post Menopause') + 
  labs(y = 'OS')

pp = pre + post  
ggsave('results/metabric/mutations/gata3_mut_type_pre_post_hr.png', width = 16, height = 8, plot = pp)
ggsave('results/metabric/mutations/gata3_mut_type_pre_post_hr.pdf', width = 16, height = 8, plot = pp)
# pdf('results/metabric/mutations/gata3_mut_type_pre_post_hr.pdf', width = 10, height = 6)
# pre
# post

os_analysis_gata3_mut_type[[1]]$Pre$plot
os_analysis_gata3_mut_type[[1]]$Post$plot
dev.off()
# ggsave(plot = pp, filename = 'results/gata3_mut_type_pre_post.pdf', width = 10, height = 6)


# CNA analysis --------
source('/orfeo/cephfs/scratch/area/vgazziero/CDSlab/brca_cro_prj/BRCA_CRO/new_analysis/analysis/CNA/plot_utils.R')

brca_drivers_genes = readRDS("../analysis/r_obj/brca_drivers.rds")
cna = readRDS('data/r_obj/metabric/data_muts/cna_selected.rds')

cna_status = cna %>% 
  full_join(. ,all_metadata, by = join_by('sample_id' == 'PATIENT_ID'))

cna_status = cna_status %>% 
  filter(!is.na(Hugo_Symbol)) %>% 
  filter(!is.na(cna_status)) %>% 
  mutate(cna_lv = case_when(
    sign(cna_status) == -1 ~ 'Deletion', 
    sign(cna_status) == 1 ~ 'Gain', 
    .default = 'None'
  )) 

cna_status_tabs = lapply(c('Gain', 'Deletion'), function(x) {
  create_table(cna_status, x)
})
names(cna_status_tabs) = c('Gain', 'Deletion')

test_results = lapply(cna_status_tabs, function(cna_type) {
  lapply(cna_type, function(x) {
    x %>%
      as.data.frame.matrix %>% 
      janitor::fisher.test()
    })
})

res_by_gene = lapply(names(test_results), function(cna_type) {
  lapply(test_results[[cna_type]] %>% names, function(g) {
    broom::tidy(test_results[[cna_type]][[g]]) %>% 
      mutate(gene = g, 
             cna_alteration = cna_type)
  }) %>% bind_rows()
}) %>% bind_rows()

saveRDS(res_by_gene, 'results/metabric/cnas/fisher_by_gene_cna_type.rds')




# old version
n_samples_tot = cna_status %>% filter(!is.na(Hugo_Symbol)) %>%filter(!is.na(cna_status)) %>% pull(sample_id) %>% unique() %>% length()

fc_by_gene = lapply(cna_status_tabs %>% names(), function(x) {
  
  print(x)
  gene = x
  
  x = cna_status_tabs[[x]] %>% 
    as.data.frame.matrix %>% 
    tibble::rownames_to_column('Menopausal_status')
  
  tot_alt = sum(x$Gain, x$Deletion)/n_samples_tot 
  
  if(tot_alt < 0.05) {
    
    res = tibble(
      gene = gene, 
      FC_gain = NA, 
      FC_del = NA
    )
    
  } else {
    
    props = x %>% 
      rowwise() %>% 
      mutate(tot = sum(None, Gain, Deletion)) %>% 
      mutate(prop_gain = Gain/tot, prop_del = Deletion/tot) %>% 
      ungroup()
    
    # FC is defined as PRE/POST
    res = props %>% 
      dplyr::select(Menopausal_status, prop_gain, prop_del) %>% 
      pivot_wider(names_from = Menopausal_status, values_from = c(prop_gain, prop_del)) %>% 
      mutate(FC_gain = log(prop_gain_Pre/prop_gain_Post, base = 2), FC_del = log(prop_del_Pre/prop_del_Post, base = 2)) %>% 
      mutate(gene = gene) %>% 
      select(gene, starts_with('FC'))
    
  }
  
  return(res)
  
}) %>% 
  bind_rows()
saveRDS(fc_by_gene, 'results/metabric/cnas/fc_by_gene_cna.rds')


test_res_all = lapply(test_results %>% names, function(x) {
  test_results[[x]] %>% 
    broom::tidy() %>% 
    mutate(gene = x) %>% 
    select(gene, everything())
}) %>% bind_rows()

x = full_join(fc_by_gene,test_res_all, by = "gene")
saveRDS(x, 'results/metabric/cnas/test_alterations_by_genes.rds')

x = readRDS('results/metabric/cnas/test_alterations_by_genes.rds')

x = x  %>% 
  mutate(is_driver = ifelse(gene %in% brca_drivers_genes, TRUE, FALSE)) %>% 
  mutate(significant = ifelse(p.value <= 0.05, 'Significant', 'ns'))
write.table(x, 'new_res/cna_alterations_by_gene_testing.csv', sep = ',', quote = F, col.names = T, row.names = F)

# x %>% 
#   filter(p.value < .05) %>% 
#   filter(gene %in% brca_drivers_genes)

del = plot_volcano_cna_genes(x, 
                       which = 'del', 
                       pth = .05, 
                       driver_list = brca_drivers_genes) + 
  xlim(c(-3.5,3.5)) + 
  ggtitle('Enrichment of deletions')
ggsave('results/metabric/cnas/enrichment_by_gene_deletions_green.png', width = 8, height = 8)
ggsave('results/metabric/cnas/enrichment_by_gene_deletions_green.pdf', width = 8, height = 8)

gain = plot_volcano_cna_genes(x, which = 'gain', pth = .05, driver_list = brca_drivers_genes) + 
  xlim(c(-3.5,3.5)) +
  ggtitle('') + 
  ggtitle('Enrichment of gains')
ggsave('results/metabric/cnas/enrichment_by_gene_gains_green.png', width = 8, height =8)
ggsave('results/metabric/cnas/enrichment_by_gene_gains_green.pdf', width = 8, height =8)

del + gain 

# checking the genes alterations

cna_status = cna %>% 
  full_join(. ,all_metadata, by = join_by('sample_id' == 'PATIENT_ID'))

cna_status = cna_status %>% 
  filter(!is.na(Hugo_Symbol)) %>% 
  # filter(!is.na(cna_status)) %>% 
  mutate(cna_lv = case_when(
    sign(cna_status) == -1 ~ 'Deletion', 
    sign(cna_status) == 1 ~ 'Gain', 
    .default = 'None'
  )) 

pre = menopause_status %>% filter(Inferred_menopause_status == 'Pre') %>% pull(sample_id) %>% unique %>% length
post = menopause_status %>% filter(Inferred_menopause_status == 'Post') %>% pull(sample_id) %>% unique %>% length

cna_cohort = cna_status %>% 
  dplyr::select(Hugo_Symbol, sample_id, Inferred_menopause_status, cna_lv) %>% 
  group_by(Inferred_menopause_status, Hugo_Symbol, cna_lv) %>% 
  mutate(n_samples = n()) %>% 
  mutate(
    prop = ifelse(
      Inferred_menopause_status == 'Pre', 
      n_samples/pre, 
      n_samples/post
    )
  )

cna_cohort %>% 
  filter(Hugo_Symbol == 'MED23') %>% 
  # filter(cna_lv == 'Gain') %>% 
  dplyr::select(Inferred_menopause_status, cna_lv, prop,n_samples) %>% 
  distinct()


cna_cohort %>% 
  filter(Hugo_Symbol %in% brca_drivers_genes) %>% 
  ggplot(aes(
    x = Hugo_Symbol, 
    fill = cna_lv
  )) + 
  geom_bar(stat = 'count', position = 'dodge') + 
  facet_wrap(~Inferred_menopause_status, scales = 'free') + 
  theme_bw() + 
  coord_flip()
   

x %>% 
  filter(gene == 'MED23')


