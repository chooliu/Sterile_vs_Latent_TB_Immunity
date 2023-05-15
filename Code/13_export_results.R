# ==============================================================================
# 13_export_results.R
# finally, summarize modeling results in Excel output
# ==============================================================================




# summarize maker inclusion / exclusion ========================================

check_inclusion <-
  function(p_in, ct_in, m_in, resultstable) {
    if (length(resultstable) == 3) {
      resultstable <- resultstable$Results
    }

    resultstable %>%
      ungroup() %>%
      filter(Panel == p_in & CellType == ct_in & Marker == m_in) %>%
      nrow() %>%
      `!=`(., 0) %>%
      if_else(., "yes", "")
  }

resultstable_marker_summary <-
  deltas_final %>%
  group_by(Panel, CellType, Marker) %>%
  summarize(
    `Included, Unstim L- v L+` = check_inclusion(Panel, CellType, Marker, resultstable_unstim_lneg_lpos),
    `Included, Unstim L- v BCG` = check_inclusion(Panel, CellType, Marker, resultstable_unstim_lneg_bcg),
    `Included, Unstim L+ v BCG` = check_inclusion(Panel, CellType, Marker, resultstable_unstim_lpos_bcg),
    `Included, Stim L- v L+` = check_inclusion(Panel, CellType, Marker, resultstable_stim_lneg_lpos),
    `Included, Stim L- v BCG` = check_inclusion(Panel, CellType, Marker, resultstable_stim_lneg_bcg),
    `Included, Stim L+ v BCG` = check_inclusion(Panel, CellType, Marker, resultstable_stim_lpos_bcg),
  )

resultstable_marker_summary %>%
  dplyr::select(contains("Included")) %>%
  ungroup() %>%
  summarize_all(function(x) {
    sum(x == "yes")
  })


resultstable_marker_summary_2 <-
  deltas_final %>%
  group_by(Panel, CellType, Marker, Group) %>%
  summarize(
    `Unstim%Zeros` = sum(Percent_NoStim == 0) / n(),
    `Unstim%LessOne` = sum(Percent_NoStim < 1) / n(),
    `UnstimMin` = min(Percent_NoStim),
    `UnstimQ25` = quantile(Percent_NoStim, 0.25),
    `UnstimMedian` = quantile(Percent_NoStim, 0.5),
    `UnstimQ75` = quantile(Percent_NoStim, 0.75),
    `UnstimMax` = max(Percent_NoStim),
    `Stim%Zeroes` = sum(Percent_CMStim == 0) / n(),
    `Stim%LessOne` = sum(Percent_CMStim < 1) / n(),
    `StimMin` = min(Percent_CMStim),
    `StimQ25` = quantile(Percent_CMStim, 0.25),
    `StimMedian` = quantile(Percent_CMStim, 0.5),
    `StimQ75` = quantile(Percent_CMStim, 0.75),
    `StimMax` = max(Percent_CMStim),
    `DiffMin` = min(Difference),
    `DiffQ25` = quantile(Difference, 0.25),
    `DiffMedian` = quantile(Difference, 0.5),
    `DiffQ75` = quantile(Difference, 0.75),
    `DiffMax` = max(Difference),
  ) %>%
  pivot_longer(cols = 5:ncol(.)) %>%
  pivot_wider(
    id_cols = c("Panel", "CellType", "Marker"),
    names_from = c("name", "Group"), names_sep = " "
  )




# export betareg results as Excel file =========================================

list(
  `Converting to Markers` = parsed_features_freq,
  `Marker Freq Summary` = bind_cols(
    resultstable_marker_summary,
    resultstable_marker_summary_2[, -c(1:3)]
  ),
  `Marker Count Summary` = parent_count_summary,
  `Unstim L- v L+` = resultstable_unstim_lneg_lpos$Results,
  `Unstim L- v BCG` = resultstable_unstim_lneg_bcg$Results,
  `Unstim L+ v BCG` = resultstable_unstim_lpos_bcg$Results,
  `Stim L- v L+` = resultstable_stim_lneg_lpos$Results,
  `Stim L- v BCG` = resultstable_stim_lneg_bcg$Results,
  `Stim L+ v BCG` = resultstable_stim_lpos_bcg$Results,
  `Stim L- v L+ (nbl)` = resultstable_stim_lneg_lpos_nbl$Results,
  `Stim L- v BCG (nbl)` = resultstable_stim_lneg_bcg_nbl$Results,
  `Stim L+ v BCG (nbl)` = resultstable_stim_lpos_bcg_nbl$Results,
  `Parent Unstim L- v L+` = resultstable_unstim_lneg_lpos_parent$Results,
  `Parent Unstim L- v BCG` = resultstable_unstim_lneg_bcg_parent$Results,
  `Parent Unstim L+ v BCG` = resultstable_unstim_lpos_bcg_parent$Results,
  `Phenog Unstim L- v L+` = resultstable_unstim_lneg_lpos_phenograph$Results %>% select(-CellType),
  `Phenog Stim L- v L+` = resultstable_stim_lneg_lpos_phenograph$Results %>% select(-CellType),
  `Phenog Stim L- v L+ (nbl)` = resultstable_stim_lneg_lpos_phenograph_nbl$Results %>% select(-CellType),
  `Phenog Manual Cluster` = results_phenograph_TCluster26
) %>%
  set_names(., paste0(1:length(.), " ", names(.))) %>%
  write_xlsx(
    path = paste0("./Reports/", current_date_label, "_TB_R21_CompiledResults.xlsx"),
    format_headers = F
  )

system(
  paste0("open ./Reports/", current_date_label, "_TB_R21_CompiledResults.xlsx"))



