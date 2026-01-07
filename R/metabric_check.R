# analysing metabric data to integrate them with the cohort analysis from cro
rm(list=ls())
library(tidyverse)
library(readxl)

# get the sequenced cases
path = 'data/brca_metabric/case_lists/'
list.files(path)

all_cases = read.table('data/brca_metabric/case_lists/cases_sequenced.txt', sep = '\t', skip = 5)
all_cases = all_cases['1', ] %>% as.character()
all_cases = gsub('case_list_ids: ', '', all_cases)
length(all_cases)
# read data from the selected metadata
metadata_selected = readxl::read_excel('data/METABRIC FILTERED.xlsx')
metadata_selected = metadata_selected %>% 
  rename(sample_id = 'Sample ID') %>% 
  rename(patient_id = 'Patient ID')
samples = metadata_selected$sample_id %>% unique

samples_sequenced = intersect(all_cases, samples)
patients = metadata_selected$patient_id %>% unique
length(patients)
length(samples_sequenced)

# filtering metadata
metadata_patient = read.table('data/brca_metabric/data_clinical_patient.txt', 
                              sep = '\t', header = T) %>% 
  as_tibble()

metadata_patient = metadata_patient %>% 
  dplyr::filter(PATIENT_ID %in% samples_sequenced)

saveRDS(metadata_selected, 'data/r_obj/metadata/metadata_selected_samples.rds')
saveRDS(metadata_patient, 'data/r_obj/metadata/metadata_selected_patients.rds')

# selecting samples and mutations
muts = read.table('data/brca_metabric/data_mutations.txt', sep = '\t', skip = 1, 
                  header = T, comment.char = "", quote = "") %>% 
  as_tibble()

muts_selected = muts %>% 
  filter(Tumor_Sample_Barcode %in% samples_sequenced) 
muts_selected$Tumor_Sample_Barcode %>% unique %>% length

intersect(muts$Tumor_Sample_Barcode, samples_sequenced) %>% length()

saveRDS(muts_selected, 'data/r_obj/metabric/data_muts/muts_selected.rds')


# selecting cnas
cna = read.table('data/brca_metabric/data_cna.txt', sep = '\t', header = T, 
                 comment.char = "", quote = "") %>% 
  as_tibble()

colnames(cna) = gsub('\\.', '-', colnames(cna))

samples_cna = intersect(samples_sequenced, colnames(cna))

cna_longer = cna %>% 
  dplyr::select(Hugo_Symbol, any_of(samples_cna)) %>% 
  pivot_longer(names_to = 'sample_id', cols = samples_cna, values_to = 'cna_status')

saveRDS(cna_longer, 'data/r_obj/metabric/data_muts/cna_selected.rds')
