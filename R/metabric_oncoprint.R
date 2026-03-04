rm(list=ls())
setwd('/orfeo/cephfs/scratch/area/vgazziero/CDSlab/brca_cro_prj/brca_extended_analysis')
library(tidyverse)
library(ggsurvfit)
library(survival)
library(wesanderson)
library(survminer)
library(ComplexHeatmap)
library(rcartocolor)
source('R/utils.R')

muts_selected = readRDS('data/r_obj/metabric/data_muts/muts_selected.rds')
brca_drivers = readRDS('../analysis/r_obj/brca_drivers.rds')

muts_selected_drivers = muts_selected %>% 
  filter(Hugo_Symbol %in% brca_drivers) %>% 
  select(Hugo_Symbol, HGVSp, Chromosome, Start_Position, End_Position, starts_with('HGVS'), Consequence, Variant_Classification) %>% 
  distinct() 

all_metadata = readRDS('data/r_obj/metabric/metadata/all_metadata.rds')

muts_selected = muts_selected %>% 
  mutate(is_driver = ifelse(Hugo_Symbol %in% brca_drivers, TRUE, FALSE)) %>% 
  mutate(is_driver = ifelse(Variant_Classification %in% c('Silent', 'Intron'), FALSE, is_driver))

# add metadata of menopause
menopause_status = all_metadata %>% 
  select(sample_id, Inferred_menopause_status)

muts_selected = muts_selected %>% 
  full_join(., menopause_status, by = join_by('Tumor_Sample_Barcode' == 'sample_id')) %>% 
  mutate(Inferred_menopause_status = factor(Inferred_menopause_status, levels = c('Pre', 'Post', 'Unknown')))

muts_selected %>% 
  filter(Hugo_Symbol == 'GATA3', is_driver == TRUE) %>% 
  group_by(Inferred_menopause_status) %>% 
  count()

cna = readRDS('data/r_obj/metabric/data_muts/cna_selected.rds')

cna_status = cna %>% 
  filter(!is.na(Hugo_Symbol)) %>% 
  filter(!is.na(cna_status)) %>% 
  mutate(cna_lv = case_when(
    sign(cna_status) == -1 ~ 'Deletion', 
    sign(cna_status) == 1 ~ 'Gain', 
    .default = 'None'
  ))

muts = muts_selected %>% 
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification, is_driver)

mut_cna = cna_status %>% 
  full_join(., muts, by = join_by(
    'sample_id' == 'Tumor_Sample_Barcode', 
    'Hugo_Symbol' == 'Hugo_Symbol'
  ))

mut_cna = mut_cna %>% 
  mutate(Variant_Classification = ifelse(is.na(Variant_Classification), 'Wild-type', Variant_Classification)) %>% 
  mutate(cna_lv = ifelse(is.na(cna_lv), 'Not estimated', cna_lv))

mut_cna_driver = mut_cna %>% 
  filter(is_driver == TRUE)

mut_cna_driver = mut_cna_driver %>% 
  group_by(sample_id, Hugo_Symbol) %>% 
  mutate(Variant_Classification = ifelse(n() > 1, 'multihit', Variant_Classification)) %>% 
  distinct() %>% 
  # mutate(mut_cna = paste(Variant_Classification, cna_lv, sep = ',')) %>% 
  mutate(mut_cna = Variant_Classification) %>% 
  select(Hugo_Symbol, sample_id, mut_cna) %>% 
  group_by(Hugo_Symbol) %>% 
  filter(n() > 30) %>% 
  pivot_wider(names_from = sample_id, values_from = mut_cna) %>% 
  tibble::column_to_rownames('Hugo_Symbol')


# define the graphics
# define the alteration function
cna_types = setNames(
  nm = c('None', 'Gain', 'Deletion', 'Not estimated'), 
  object = c(8, 17, 18, NA)
)

cn_fun = lapply(cna_types, function(s) {
  function(x, y, w, h) 
    grid.points(x,y,w,h, pch = unname(s), size = unit(3, "mm"))
})

mut_types = c("Missense_Mutation",
              "Nonsense_Mutation",
              "multihit",
              "Frame_Shift_Del",
              "Splice_Site",
              "Frame_Shift_Ins",
              "In_Frame_Ins",
              "Splice_Region",
              "In_Frame_Del",
              "Nonstop_Mutation",
              "Translation_Start_Site")

mut_cols = setNames(nm = mut_types, 
         object = adjustcolor(carto_pal(length(c(mut_types, "NA")), "Prism"), alpha.f = 0.9)[1:length(mut_types)])

vars_fun = lapply(mut_cols, function(x) {
  alter_graphic("rect", width = 0.8, height = 0.8, fill = unname(x))
})
vars_fun = c(vars_fun, 'background' = alter_graphic("rect", width = 0.8, height = 0.8, fill = "#CCCCCC"))

alt_fun = c(vars_fun, cn_fun)

# annotation 

top_annotation_colors = list(
  "Inferred_menopause_status" = setNames(
    nm = c('Pre', 'Post'),
    # object = c("goldenrod1", "#CBA1D2")
    object = c(rgb(210, 231, 167, maxColorValue = 255), '#CBA1D2')
  ))

meno_status = menopause_status %>% 
  filter(sample_id %in% colnames(mut_cna_driver)) %>% 
  select(Inferred_menopause_status)

bottom_annotation = create_annotation(x = meno_status, 
                                      ann_colors = top_annotation_colors, 
                                      position = "column")

# might remove the cnas?
png('results/metabric/mutations/oncoprint_drivers.png', width = 20, height = 12, units = 'cm', res = 300)
oc = oncoPrint(mut_cna_driver %>% as.matrix(), 
          alter_fun = alt_fun, 
          remove_empty_columns = T, 
          remove_empty_rows = T,
          show_row_names = T, 
          col = mut_cols, 
          pct_gp = gpar(fontsize = 7),
          row_names_gp = gpar(fontsize = 7), 
          name = 'Alteration types',
          top_annotation = bottom_annotation,
          column_split = meno_status %>% pull(Inferred_menopause_status),
          heatmap_legend_param = list(direction = "horizontal", nrow = 3))
draw(oc, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()
