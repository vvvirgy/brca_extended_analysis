rm(list=ls())
.libPaths(new = '/u/area/vgazziero/R/rstudio_4_4/')
library(tidyverse)
library(janitor)

# CNA analysis --------
source('/orfeo/cephfs/scratch/area/vgazziero/CDSlab/brca_cro_prj/BRCA_CRO/new_analysis/analysis/CNA/plot_utils.R')
source('/orfeo/cephfs/scratch/area/vgazziero/CDSlab/brca_cro_prj/brca_extended_analysis/R/utils.R')

brca_drivers_genes = readRDS("../analysis/r_obj/brca_drivers.rds")
cna = readRDS('data/r_obj/metabric/data_muts/cna_selected.rds')
all_metadata = readRDS('data/r_obj/metabric/metadata/all_metadata.rds')

# add metadata of menopause
menopause_status = all_metadata %>% 
  dplyr::select(sample_id, Inferred_menopause_status)

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

res_by_gene = res_by_gene %>% 
  mutate(estimate = log(estimate, 2)) %>%
  rename(log2FC = estimate) %>%
  mutate(is_driver = ifelse(gene %in% brca_drivers_genes, TRUE, FALSE)) %>% 
  mutate(significant = ifelse(p.value <= 0.05, 'Significant', 'ns'))

plot_volcano_cna_genes(res_by_gene, which = 'Deletion', driver_list = brca_drivers_genes)
plot_volcano_cna_genes(res_by_gene, which = 'Gain', driver_list = brca_drivers_genes)

plot_volcano_cna_genes(res_by_gene, 
                       which = c('Deletion', 'Gain'), 
                       driver_list = brca_drivers_genes, 
                       cols = menopause_cols) + 
  # xlim(c(-3.5,3.5)) + 
  ggtitle('Enrichment of CNA alterations')

ggsave('results/metabric/cnas/enrichment_by_gene.png', width = 12, height = 8)
ggsave('results/metabric/cnas/enrichment_by_gene.pdf', width = 12, height = 8)



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
  filter(Hugo_Symbol == 'NOTCH1') %>% 
  # filter(cna_lv == 'Gain') %>% 
  dplyr::select(Inferred_menopause_status, cna_lv, prop,n_samples) %>% 
  distinct() %>% 
  filter(cna_lv == 'Deletion')


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


