# test for enrichment of cna events between different groups

create_table = function(data, 
                        test_lv) {
  data %>% 
    mutate(cna_lv = ifelse(cna_lv != test_lv, 'Other', test_lv)) %>% 
    dplyr::count(Hugo_Symbol, Inferred_menopause_status, cna_lv) %>%
    complete(Hugo_Symbol, Inferred_menopause_status, cna_lv, fill = list(n = 0)) %>% 
    group_by(Hugo_Symbol) %>% 
    summarise(
      tab = list(xtabs(n ~ Inferred_menopause_status + cna_lv, data = cur_data(), drop.unused.levels = FALSE))
    ) %>% 
    deframe()
  
}
