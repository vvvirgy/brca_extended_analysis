
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
    group_by(Inferred_menopause_status, Hugo_Symbol) %>% 
    mutate(n = n(), prop = n/tot_samples) %>% 
    select(Hugo_Symbol, Inferred_menopause_status, n, tot_samples, prop) %>% 
    distinct()
}

plot_prop_samples = function(df, colors) {
  df %>% 
    ggplot(aes(y = Hugo_Symbol, fill = Inferred_menopause_status, x = prop)) + 
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
    ggplot(aes(y = Hugo_Symbol, fill = Inferred_menopause_status, x = prop)) + 
    geom_bar(stat = 'identity') + 
    theme_bw() + 
    scale_x_continuous(labels = scales::percent) +
    scale_fill_manual(values = colors) + 
    labs(x = 'Percentage of samples mutated', y = 'Gene symbols') + 
    theme(legend.position = 'right') + 
    guides(fill = guide_legend(title = 'Menopause status', nrow = 3))
}


# survival analysis ------- 
survival_analysis = function(data, time, event, what, gene) {
  
  data = data %>% 
    dplyr::filter(!is.na(!!gene))
  
  surv_obj <- Surv(time = (data %>% dplyr::pull(time)), event = (data %>% dplyr::pull(event)) )
  
  # Create formula dynamically
  # formula_obj <- as.formula(paste("surv_obj ~", what, '+ Inferred_menopause_status'))
  # formula_obj2 <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", what, '+ Inferred_menopause_status'))
  formula_obj <- as.formula(paste("surv_obj ~", what))
  formula_obj2 <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", what))
  
  # Fit
  fit_km <- survfit(formula_obj, data = data)
  diff <- survdiff(formula_obj, data = data)
  
  # perform also cox fit to have the comparison among each category 
  # data_coxph <- data %>%
  #   mutate(surv_object = Surv(  time  = .data[[time]],
  #                               event = .data[[event]]))
  
  data_df = as.data.frame(data)

  fit_coxph = coxph(formula = formula_obj2, ties = 'exact', data = data_df)
  
  hr_plot = tryCatch(expr = {
    ggforest(fit_coxph, data = data_df)}, 
    error = function(e) {
      ggplot()
    })
  
  
  # summary(fit_coxph)
  # coxfit <- coxph(
  #   Surv(os_years, os_status) ~ what,
  #   data = data,
  #   ties = 'exact')
  
  
  # plot
  pp = survfit2(formula_obj2, data = data) %>% 
    ggsurvfit(linewidth = 1.5) +
    labs( x = 'Time (years)', y = 'probs')  + 
    # add_confidence_interval() +
    scale_color_brewer(palette = "Set1") +
    add_quantile() + 
    # scale_color_manual(values = wes_palette(3, name = "Darjeeling1", type = "discrete")) +
    # scale_fill_manual(values = wes_palette(3, name = "Darjeeling1", type = "discrete")) + 
    add_pvalue(location = 'annotation') +
    ggtitle(what) + 
    theme(axis.text = element_text(family = "Arial", size = 13),
          axis.title = element_text(family = "Arial", size = 13),
          legend.text = element_text(family = "Arial", size = 12),
          panel.border = element_rect(linewidth = 1.5), 
          title = element_text(size = 15), 
          legend.direction = 'vertical') + 
    guides(color = guide_legend(ncol = 2))
  
  list('fit' = fit_km, 
       'Diff' = diff, 
       'coxph' = fit_coxph,
       'hr_plot' = hr_plot,
       'km_plot' = pp)
}


# oncoprint ----
create_annotation = function(x, ann_colors, position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position, 
                                    show_annotation_name = F, 
                                    annotation_legend_param = list(nrow = 2, ncol = 3, width = 12, by_row = T))
  
}

# cna analysis ----

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

# plot utils 

plot_volcano_cna_genes = function(x, which, pth = .05, driver_list = NULL, cols) {
  
  df = x %>% 
    dplyr::select(gene, log2FC, p.value, cna_alteration) %>% 
    # filter(if_any(starts_with("FC"), ~ !is.infinite(.x)))
    filter(!is.na(log2FC)) %>% 
    filter(!is.infinite(log2FC))
  
  # add classes
  df = df %>% 
    mutate(cls = 
             case_when(
               p.value <= pth ~ 'Significant', 
               .default = 'ns'
             ))  %>% 
    mutate(log2FC = log2FC*-1) %>% 
    filter(cna_alteration %in% which)
  
  # add driver labels
  df = df %>% 
    mutate(is_driver = ifelse(gene %in% driver_list, TRUE, FALSE)) %>% 
    mutate(is_driver = factor(is_driver, levels = c(TRUE, FALSE))) %>% 
    mutate(dr_label = ifelse((is_driver == TRUE & cls == 'Significant'), gene, NA))
  
  df %>% 
    # filter(is_driver) %>% 
    ggplot(aes(x = log2FC, 
               y = -log(p.value, base = 10), 
               shape = is_driver, 
               label = dr_label,
               # size = is_driver, 
               color = cls)) + 
    geom_point() +
    # geom_label() +
    scale_color_manual(values = c('Significant' = 'navyblue', 'ns' = 'grey20')) +
    scale_shape_manual(values = setNames(object = c(8, 1), c(TRUE, FALSE))) + 
    geom_hline(yintercept = -log(pth, base = 10), 
               linetype = 'dashed', 
               colour = 'firebrick'
    ) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed', 
               colour = 'grey'
    ) + 
    theme_bw() + 
    ggplot2::geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf, Sample_type = 'post'),
                       aes(xmin = xmin, xmax =xmax, ymin = ymin, ymax = ymax, fill = Sample_type),
                       alpha=0.3, inherit.aes = F)+
    ggplot2::geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf, Sample_type = 'pre'),
                       aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = Sample_type),
                       alpha=0.3, inherit.aes = F) + 
    # ggplot2::scale_fill_manual(values = c('pre' = "goldenrod1", 'post' = '#CBA1D2')) +
    ggplot2::scale_fill_manual(values = c('pre' = unname(cols['Pre']), 'post' = unname(cols['Post'])), 
      # values = c('pre' = rgb(210, 231, 167, maxColorValue = 255), 'post' = '#CBA1D2'),
                               labels = c('pre' = 'Pre−Menopause', 'post' = 'Post−Menopause')) +
    ggrepel::geom_text_repel(max.overlaps = 100, show.legend = F, color = 'black') +
    ylab('-log(pVal)') + 
    facet_wrap(~cna_alteration, scales = 'free')
}

# plot volcano mutations
plot_volcano_mut_genes = function(x, pth = .05, driver_list = NULL, cols) {
  
  df = x %>% 
    dplyr::select(gene, log2FC, p.value) %>% 
    # filter(if_any(starts_with("FC"), ~ !is.infinite(.x)))
    filter(!is.na(log2FC)) %>% 
    filter(!is.infinite(log2FC))
  
  # add classes
  df = df %>% 
    mutate(cls = 
             case_when(
               p.value <= pth ~ 'Significant', 
               .default = 'ns'
             ))  %>% 
    mutate(log2FC = log2FC*-1) 
  
  # add driver labels
  df = df %>% 
    mutate(is_driver = ifelse(gene %in% driver_list, TRUE, FALSE)) %>% 
    mutate(is_driver = factor(is_driver, levels = c(TRUE, FALSE))) %>% 
    mutate(dr_label = ifelse((is_driver == TRUE & cls == 'Significant'), gene, NA))
  
  df %>% 
    # filter(is_driver) %>% 
    ggplot(aes(x = log2FC, 
               y = -log(p.value, base = 10), 
               shape = is_driver, 
               label = dr_label,
               # size = is_driver, 
               color = cls)) + 
    geom_point() +
    # geom_label() +
    scale_color_manual(values = c('Significant' = 'navyblue', 'ns' = 'grey20')) +
    scale_shape_manual(values = setNames(object = c(8, 1), c(TRUE, FALSE))) + 
    geom_hline(yintercept = -log(pth, base = 10), 
               linetype = 'dashed', 
               colour = 'firebrick'
    ) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed', 
               colour = 'grey'
    ) + 
    theme_bw() + 
    ggplot2::geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf, Sample_type = 'post'),
                       aes(xmin = xmin, xmax =xmax, ymin = ymin, ymax = ymax, fill = Sample_type),
                       alpha=0.3, inherit.aes = F)+
    ggplot2::geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf, Sample_type = 'pre'),
                       aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = Sample_type),
                       alpha=0.3, inherit.aes = F) + 
    # ggplot2::scale_fill_manual(values = c('pre' = "goldenrod1", 'post' = '#CBA1D2')) +
    ggplot2::scale_fill_manual(values = c('pre' = unname(cols['Pre']), 'post' = unname(cols['Post'])), 
                               # values = c('pre' = rgb(210, 231, 167, maxColorValue = 255), 'post' = '#CBA1D2'),
                               labels = c('pre' = 'Pre−Menopause', 'post' = 'Post−Menopause')) +
    ggrepel::geom_text_repel(max.overlaps = 100, show.legend = F, color = 'black') +
    ylab('-log(pVal)') 
}


