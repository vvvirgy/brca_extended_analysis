
# Colors plotting ------ 

# colors menopause
menopause_cols = setNames(
  nm = c('Pre', 'Post', 'Unknown'),
  # object = c("goldenrod1", "#CBA1D2")
  # object = c(rgb(210, 231, 167, maxColorValue = 255), '#CBA1D2')
  object = c('#ffc125', '#CBA1D2', '#EBE1D1')
)

# utils functions ----- 

# create proportions

get_props = function(df, group) {
  df %>% 
    group_by(across(all_of(group))) %>% 
    mutate(tot_samples = n()) %>%
    group_by(INFERRED_MENOPAUSAL_STATE, Hugo_Symbol) %>% 
    mutate(n = n(), prop = n/tot_samples) %>% 
    select(Hugo_Symbol, INFERRED_MENOPAUSAL_STATE, n, tot_samples, prop) %>% 
    distinct()
}

plot_prop_samples = function(df, colors) {
  df %>% 
    ggplot(aes(y = Hugo_Symbol, fill = INFERRED_MENOPAUSAL_STATE, x = prop)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    theme_bw() + 
    scale_x_continuous(labels = scales::percent) + 
    scale_fill_manual(values = colors) + 
    labs(x = 'Percentage of samples', y = 'Gene symbols') + 
    theme(legend.position = 'right') + 
    guides(fill = guide_legend(title = 'Menopause status', nrow = 3))
}

plot_prop_genes = function(df, colors) {
  df %>% 
    ggplot(aes(y = Hugo_Symbol, fill = INFERRED_MENOPAUSAL_STATE, x = prop)) + 
    geom_bar(stat = 'identity') + 
    theme_bw() + 
    scale_x_continuous(labels = scales::percent) +
    scale_fill_manual(values = colors) + 
    labs(x = 'Percentage of samples mutated', y = 'Gene symbols') + 
    theme(legend.position = 'right') + 
    guides(fill = guide_legend(title = 'Menopause status', nrow = 3))
}
