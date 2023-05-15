# ==============================================================================
# 12_targeted_marker_viz.R
# generic functions for plotting all % positivity by group
# ==============================================================================




# format LTBI- vs + residuals ==================================================
# in long format for plotting

predicted_table_to_long <-
  function(resultsobj, dat_start_index = 5) {
    resultsobj$Residuals %>%
      pivot_longer(
        cols = all_of(dat_start_index):ncol(.),
        names_to = "PID", values_to = "OutcomeOfInterest"
      ) %>%
      left_join(.,
        deltas_final %>% filter(!duplicated(PID)) %>%
          select(PID, Group, Run),
        by = "PID"
      )
  }

predictedvals_unstim_long <-
  predicted_table_to_long(resultstable_unstim_lneg_lpos)
predictedvals_stim_long <-
  predicted_table_to_long(resultstable_stim_lneg_lpos)
predictedvals_stim_nbl_long <-
  predicted_table_to_long(resultstable_stim_lneg_lpos_nbl)




# generic marker plotting fxn ==================================================
# plots residualized values by default

plot_activation_marker <-
  function(panel_in, celltype_in, marker_in,
           inputtable = NULL, ylabel = NULL,
           ysuffix = "\n (Residualized % Positivity, Unstim)",
           export_prefix = "",
           exclude_BCG = F,
           add_hline = F) {
    
    if (is.null(ylabel)) {
      ylabel <-
        paste(c(panel_in, celltype_in, marker_in), collapse = " // ") %>%
        paste(ysuffix)
    }

    xlabels <- c("BCG", "Resist.", "LTBI")

    inputtable <-
      inputtable %>%
      filter(Panel == panel_in & CellType == celltype_in & Marker == marker_in)

    if (exclude_BCG) {
      inputtable <-
        inputtable %>%
        filter(Group != "BCG")
      xlabels <- xlabels[-1]
    }

    plot_out <-
      ggplot(
        inputtable,
        aes(x = Group, y = OutcomeOfInterest, color = Group)
      ) +
      geom_quasirandom(alpha = 0.6, size = 2, width = 0.2) +
      stat_summary(fun.data = fxnMakeQuantiles,
                   geom = "errorbar", color = "black",
                   width = 0.3, size = 0.8, alpha = 0.8) +
      stat_summary(fun = median, geom = "point", color = "black",
                   shape = 15, size = 3) +
      theme_few(base_size = 10) +
      scale_y_continuous(ylabel, breaks = pretty_breaks(n = 6)) +
      scale_color_manual(values = palette_cohorts) +
      theme(legend.position = "none") +
      scale_x_discrete(labels = xlabels)

    if (add_hline) {
      plot_out <-
        plot_out +
        geom_hline(yintercept = 0, lty = 3)
    }

    if (export_prefix != "") {
      ggsave(
        plot = plot_out,
        filename = paste0(
          "./FiguresTables/", current_date_label, "/",
          export_prefix,
          panel_in, "_", celltype_in, "_", marker_in, ".png"
        ),
        width = 3, height = 3
      )
    }

    plot_out
  }




# plotting examples, markers ===================================================

# single marker example
plot_activation_marker(
  "APC", "cDC1", "IL10",
  predictedvals_unstim_long
)

# example of checking raw values, incl. for run-specific effects
plot_activation_marker("APC", "cDC1", "IL10",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_NoStim),
  ysuffix = "\n Raw % Positivity, Unstim",
  export_prefix = "rawunstim"
) +
  facet_grid(~Run)

# plot top from stim LTBI+/- results (baseline adjusted)
plot_activation_marker(resultstable_stim_lneg_lpos$Results$Panel[1],
  resultstable_stim_lneg_lpos$Results$CellType[1],
  resultstable_stim_lneg_lpos$Results$Marker[1],
  deltas_final %>% mutate(OutcomeOfInterest = Difference),
  ysuffix = "\n(Raw % of Events, Stim)"
) +
  facet_grid(~Run)

# or stim, without baseline adjustment ("nbl")
plot_activation_marker(resultstable_stim_lneg_lpos_nbl$Results$Panel[1],
  resultstable_stim_lneg_lpos_nbl$Results$CellType[1],
  resultstable_stim_lneg_lpos_nbl$Results$Marker[1],
  predictedvals_stim_nbl_long,
  ysuffix = "\n(Residualized % of Events, Stim) [nbl]"
)


# all signif results requested for paper, revision 1 ===========================

# loop through all signif unstim markers
dir.create("./FiguresTables/20220925/Scatterplot_Unstim/")
resultstable_unstim_lneg_lpos$Results %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  mutate(ylabel = paste0(CellType, " ", Marker, "\n(Unstimulated %)")) %>%
  pmap(function(Panel, CellType, Marker, ylabel) {
    plot_activation_marker(Panel, CellType, Marker,
                           inputtable = deltas_final %>% mutate(OutcomeOfInterest = Percent_NoStim),
                           export_prefix = "Scatterplot_Unstim/",
                           ylabel = ylabel, exclude_BCG = T
    )
  })

# loop through all signif stim markers
dir.create("./FiguresTables/20220925/Scatterplot_Delta/")
resultstable_stim_lneg_lpos$Results %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  mutate(ylabel = paste0(CellType, " ", Marker, "\n(Δ: stim % - unstim %)")) %>%
  pmap(function(Panel, CellType, Marker, ylabel) {
    plot_activation_marker(Panel, CellType, Marker, ylabel = ylabel,
                           inputtable = deltas_final %>%
                             mutate(OutcomeOfInterest = Difference),
                           export_prefix = "Scatterplot_Delta/", exclude_BCG = T
    )
  })

# + BCG
dir.create("./FiguresTables/20220925/Scatterplot_Delta_plusBCG/")
bind_rows(
  resultstable_stim_lpos_bcg$Results,
  resultstable_stim_lneg_bcg$Results) %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  mutate(ylabel = paste0(CellType, " ", Marker, "\n(Δ: stim % - unstim %)")) %>%
  pmap(function(Panel, CellType, Marker, ylabel) {
    plot_activation_marker(Panel, CellType, Marker, ylabel = ylabel,
                           inputtable = deltas_final %>% mutate(OutcomeOfInterest = Difference),
                           export_prefix = "Scatterplot_Delta_plusBCG/", exclude_BCG = F
    )
  })

# loop through all signif stim markers (nbl)
dir.create("./FiguresTables/20220925/Scatterplot_Stim/")
resultstable_stim_lneg_lpos_nbl$Results %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  mutate(ylabel = paste0(CellType, " ", Marker, "\nMtb-stimulated %")) %>%
  pmap(function(Panel, CellType, Marker, ylabel) {
    plot_activation_marker(Panel, CellType, Marker, ylabel = ylabel,
                           inputtable = deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
                           export_prefix = "Scatterplot_Stim/", exclude_BCG = T
    )
  })

# + BCG (nbl)
dir.create("./FiguresTables/20220925/Scatterplot_Stim_plusBCG/")
bind_rows(
  resultstable_stim_lpos_bcg_nbl$Results,
  resultstable_stim_lneg_bcg_nbl$Results) %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  mutate(ylabel = paste0(CellType, " ", Marker, "\nMtb-stimulated %")) %>%
  pmap(function(Panel, CellType, Marker, ylabel) {
    plot_activation_marker(Panel, CellType, Marker, ylabel = ylabel,
                           inputtable = deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
                           export_prefix = "Scatterplot_Stim_plusBCG/", exclude_BCG = F
    )
  })



# phenograph unstim, with and without BCG group --------------------------------

plot_activation_marker("TCell", "Phenograph", "CD8+GMM+GranzB+",
  dat_deltas_phenograph_TCluster26 %>% mutate(OutcomeOfInterest = Percent_NoStim),
  export_prefix = "Targeted_Pheno_",
  ylabel = "Phenograph TCell Cluster 26\n(CD8+GMM+GranzB+; unstimulated %)",
  ysuffix = "", exclude_BCG = F
)

plot_activation_marker("TCell", "Phenograph", "CD8+GMM+GranzB+",
  dat_deltas_phenograph_TCluster26 %>% mutate(OutcomeOfInterest = Percent_NoStim),
  export_prefix = "Targeted_Pheno_nobcg_",
  ylabel = "Phenograph TCell Cluster 26\n(CD8+GMM+GranzB+; unstimulated %)",
  ysuffix = "", exclude_BCG = T
)


