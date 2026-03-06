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
all_metadata = readRDS('data/r_obj/metabric/metadata/all_metadata.rds')

# identying drivers
muts_selected = muts_selected %>% 
  mutate(is_driver = ifelse(Hugo_Symbol %in% brca_drivers, TRUE, FALSE)) %>% 
  mutate(is_driver = ifelse(Variant_Classification %in% c('Silent', 'Intron'), FALSE, is_driver))

# add metadata of menopause
menopause_status = all_metadata %>% 
  dplyr::select(sample_id, Inferred_menopause_status)

muts_selected = muts_selected %>% 
  full_join(., menopause_status, by = join_by('Tumor_Sample_Barcode' == 'sample_id')) %>% 
  mutate(Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post', 'Unknown'))) 

nPre = menopause_status %>% filter(Inferred_menopause_status == 'Pre') %>% pull(sample_id) %>% unique %>% length
nPost = menopause_status %>% filter(Inferred_menopause_status == 'Post') %>% pull(sample_id) %>% unique %>% length

mut_counts = muts_selected %>% 
  mutate(Inferred_menopause_status = as.character(Inferred_menopause_status)) %>% 
  filter(Inferred_menopause_status != 'Unkwon') %>% 
  mutate(mut_status = ifelse(
    is.na(Variant_Classification) | Variant_Classification == "Silent",
    "wild-type",
    "mutated"
  )) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Inferred_menopause_status, mut_status) %>%
  distinct() %>%
  dplyr::count(Hugo_Symbol, Inferred_menopause_status, mut_status) %>%
  filter(Inferred_menopause_status != 'Unkwon') %>% 
  tidyr::complete(
    Hugo_Symbol,
    Inferred_menopause_status,
    mut_status,
    fill = list(n = 0)
  ) %>% 
  pivot_wider(names_from = mut_status, values_from = n) %>% 
  mutate(
    `wild-type` = case_when(
      (Inferred_menopause_status == 'Pre') ~ nPre - mutated, 
      (Inferred_menopause_status == 'Post') ~ nPost - mutated, 
      .default = `wild-type`
    )
  ) %>%
  pivot_longer(cols = c(mutated, `wild-type`), names_to = 'mut_status', values_to = 'n')
  
mut_tb = mut_counts %>%
  mutate(Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post'))) %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(
    tab = list(xtabs(n ~ Inferred_menopause_status + mut_status, data = cur_data(), drop.unused.levels = FALSE))
  ) %>% 
  deframe()

mut_res = lapply(mut_tb, function(x) {
  x %>%
    as.data.frame.matrix %>% 
    janitor::fisher.test()
})

# i am looking at the odds of being mutated in the Pre wrt Post!

res_by_gene = lapply(names(mut_res), function(g) {
  
    broom::tidy(mut_res[[g]]) %>% 
      mutate(gene = g)
  }) %>% bind_rows()

saveRDS(res_by_gene, 'results/metabric/cnas/fisher_by_gene_mutations.rds')

res_by_gene = res_by_gene %>% 
  mutate(estimate = log(estimate, 2)) %>%
  rename(log2FC = estimate) %>%
  mutate(is_driver = ifelse(gene %in% brca_drivers, TRUE, FALSE)) %>% 
  mutate(significant = ifelse(p.value <= 0.05, 'Significant', 'ns'))

plot_volcano_mut_genes(res_by_gene, driver_list = brca_drivers, cols = menopause_cols)
ggsave('results/metabric/mutations/enrichment_by_gene.png', width = 8, height = 8)
ggsave('results/metabric/mutations/enrichment_by_gene.pdf', width = 8, height = 8)




