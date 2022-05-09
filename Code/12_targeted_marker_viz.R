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
          export_prefix, "_",
          panel_in, "_", celltype_in, "_", marker_in, ".png"
        ),
        width = 3, height = 3
      )
    }

    plot_out
  }




# plot significant activation markers ==========================================

# plot top from unstim LTBI+/- results
# resultstable_unstim_lneg_lpos

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

# loop through all signif unstim markers
resultstable_unstim_lneg_lpos$Results %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  pmap(function(Panel, CellType, Marker) {
    plot_activation_marker(Panel, CellType, Marker,
      inputtable = predictedvals_unstim_long,
      export_prefix = "residunstim"
    )
  })



# plot top from stim LTBI+/- results (baseline adjusted) -----------------------

# example of plotting one residualized stim marker
plot_activation_marker(resultstable_stim_lneg_lpos$Results$Panel[1],
  resultstable_stim_lneg_lpos$Results$CellType[1],
  resultstable_stim_lneg_lpos$Results$Marker[1],
  predictedvals_stim_long,
  ysuffix = "\n(Residualized % of Events, Stim)"
)

# loop through all signif stim markers
resultstable_stim_lneg_lpos$Results %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  pmap(function(Panel, CellType, Marker) {
    plot_activation_marker(Panel, CellType, Marker,
      inputtable = predictedvals_stim_long,
      export_prefix = "residstim",
      ysuffix = "\n(Residualized % of Events, Stim)"
    )
  })

# example of plotting stim (with baseline adjustment) raw data
plot_activation_marker(resultstable_stim_lneg_lpos$Results$Panel[1],
  resultstable_stim_lneg_lpos$Results$CellType[1],
  resultstable_stim_lneg_lpos$Results$Marker[1],
  deltas_final %>% mutate(OutcomeOfInterest = Difference),
  ysuffix = "\n(Raw % of Events, Stim)"
) +
  facet_grid(~Run)




# or stim, without baseline adjustment ("nbl") ---------------------------------

# example of plotting one residualized stim (nbl) marker
plot_activation_marker(resultstable_stim_lneg_lpos_nbl$Results$Panel[1],
  resultstable_stim_lneg_lpos_nbl$Results$CellType[1],
  resultstable_stim_lneg_lpos_nbl$Results$Marker[1],
  predictedvals_stim_nbl_long,
  ysuffix = "\n(Residualized % of Events, Stim) [nbl]"
)

# loop through all signif stim (nbl) markers
resultstable_stim_lneg_lpos_nbl$Results %>%
  filter(`Signif (FDR < 0.05)` == "*") %>%
  select(Panel, CellType, Marker) %>%
  pmap(function(Panel, CellType, Marker) {
    plot_activation_marker(Panel, CellType, Marker,
      inputtable = predictedvals_stim_nbl_long,
      export_prefix = "residstimNBL",
      ysuffix = "\n(Residualized % of Events, Stim) [nbl]"
    )
  })

# example of plotting stim (nbl) raw data
plot_activation_marker(resultstable_stim_lneg_lpos_nbl$Results$Panel[1],
  resultstable_stim_lneg_lpos_nbl$Results$CellType[1],
  resultstable_stim_lneg_lpos_nbl$Results$Marker[1],
  deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
  export_prefix = "rawstimNBL",
  ysuffix = "\n(Raw % of Events, Stim) [nbl]"
) +
  facet_grid(~Run)




# final requested for paper, from AW March/April 2022 ==========================

# CD4 and CD8 GranzB stim and unstim -------------------------------------------

plot_activation_marker("TCell", "Tconv CD4+", "GranzB",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_NoStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Unstim_",
  ylabel = "Tconv CD4+ GranzB+\n(unstimulated %)", ysuffix = ""
)

plot_activation_marker("TCell", "Tconv CD4+", "GranzB",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Stim_",
  ylabel = "Tconv CD4+ GranzB+\n(Mtb-stimulated %)", ysuffix = ""
)

plot_activation_marker("TCell", "Tconv CD8+", "GranzB",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_NoStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Unstim_",
  ylabel = "Tconv CD8+ GranzB+\n(unstimulated %)", ysuffix = ""
)

plot_activation_marker("TCell", "Tconv CD8+", "GranzB",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Stim_",
  ylabel = "Tconv CD8+ GranzB+\n(Mtb-stimulated %)", ysuffix = ""
)

# CD4 and CD8 PD1 stim and unstim ----------------------------------------------

plot_activation_marker("TCell", "Tconv CD4+", "PD1",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_NoStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Unstim_",
  ylabel = "Tconv CD4+ PD1+\n(unstimulated %)", ysuffix = ""
)

plot_activation_marker("TCell", "Tconv CD8+", "PD1",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_NoStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Unstim_",
  ylabel = "Tconv CD8+ PD1+\n(unstimulated %)", ysuffix = ""
)

plot_activation_marker("TCell", "Tconv CD4+", "PD1",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Stim_",
  ylabel = "Tconv CD4+ PD1+\n(Mtb-stimulated %)", ysuffix = ""
)

plot_activation_marker("TCell", "Tconv CD8+", "PD1",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Stim_",
  ylabel = "Tconv CD8+ PD1+\n(Mtb-stimulated %)", ysuffix = ""
)

# iNKT CD8 IFNg delta ----------------------------------------------------------

plot_activation_marker("TCell", "iNKT CD8+", "IFNγ",
  deltas_final %>% mutate(OutcomeOfInterest = Difference),
  exclude_BCG = T,
  export_prefix = "Targeted_Delta_",
  ylabel = "iNKT CD8+ IFNγ+\n(Δ: stim % - unstim %)",
  ysuffix = "", add_hline = T
)

# NKT CD107 stim and delta -----------------------------------------------------

plot_activation_marker("TCell", "NK", "CD107a",
  deltas_final %>% mutate(OutcomeOfInterest = Percent_CMStim),
  exclude_BCG = T,
  export_prefix = "Targeted_Stim_",
  ylabel = "NKT CD107a+\n(Mtb-stimulated %)",
  ysuffix = "", add_hline = T
)

plot_activation_marker("TCell", "NK", "CD107a",
  deltas_final %>% mutate(OutcomeOfInterest = Difference),
  exclude_BCG = T,
  export_prefix = "Targeted_Delta_",
  ylabel = "NKT CD107a+ Mtb-memory\n(Δ: stim % - unstim %)",
  ysuffix = "", add_hline = T
)

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

# delta GD CD25+ ---------------------------------------------------------------

plot_activation_marker("TCell", "γδT", "CD25",
  deltas_final %>% mutate(OutcomeOfInterest = Difference),
  exclude_BCG = T,
  export_prefix = "Targeted_Delta_",
  ylabel = "γδT CD25+\n(Δ: stim % - unstim %)",
  ysuffix = "", add_hline = T
)



