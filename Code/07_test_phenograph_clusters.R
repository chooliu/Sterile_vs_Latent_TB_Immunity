# ==============================================================================
# 07_test_phenograph_clusters.R
# load into same format as activation markers, to work a la dat_deltas
# ==============================================================================




# load phenograph cluster %s by sample =========================================
# clustering performed/exported by EJ within FlowJo, exported in .xlsx

tcell_phenograph_percents_long <-
  read_excel("DataRaw/20210709/20210707 TB Tcell merged clusters.xlsx",
             sheet = 1, skip = 1) %>%
  dplyr::select(PID, Stim, contains("cluster")) %>%
  mutate(Stim = if_else(grepl("CM", Stim), "CMStim", "NoStim")) %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Marker", values_to = "Percent") %>%
  pivot_wider(., id_cols = c("PID", "Marker"), names_prefix = "Percent_",
              names_from = Stim, values_from = Percent) %>%
  mutate(Panel = "TCell", CellType = "Phenograph") %>%
  left_join(., metadata_consolated, by = c("PID") )

apc_phenograph_percents_long <-
  read_excel("DataRaw/20210709/20210708 TB APC concat.xlsx",
             sheet = 1, skip = 1) %>%
  dplyr::select(PID, Stim, contains("cluster")) %>%
  mutate(Stim = if_else(grepl("CM", Stim), "CMStim", "NoStim")) %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Marker", values_to = "Percent") %>%
  pivot_wider(., id_cols = c("PID", "Marker"), names_prefix = "Percent_",
              names_from = Stim, values_from = Percent) %>%
  mutate(Panel = "APC", CellType = "Phenograph") %>%
  left_join(., metadata_consolated, by = c("PID") )




# join w metadata ==============================================================

deltas_phenograph <-
  bind_rows(tcell_phenograph_percents_long, apc_phenograph_percents_long) %>%
  bind_cols(FreqOf = "Parent",
            Exclude_Marker = F,
            Exclude_CellType = F) %>%
  mutate(Difference = Percent_CMStim - Percent_NoStim) %>%
  mutate(Marker = gsub("cluster |Cluster ", "", Marker) %>%
           as.numeric) %>%
  arrange(Marker) %>%
  mutate(Marker = paste0("Cluster ", Marker) %>%
           as.factor %>% fct_inorder(),
         GroupShort = fct_recode(Group,
                                 `L-` = "LTBI-", `L+` = "LTBI+", B = "BCG"))

# 64 total clusters (separate by panel)
deltas_phenograph %>%
  mutate(Cluster = paste0(.$Panel, .$Marker)) %>%
  .$Cluster %>% unique %>% length




# plot by group / run a la script 02 ===========================================

c(
  "Phenograph_Unstim", "Phenograph_Delta", "Phenograph_Stim"
) %>%
  paste0("./FiguresTables/", current_date_label, "/", .) %>%
  sapply(dir.create)

# unstim by group
deltas_phenograph %>%
  group_by(Panel) %>%
  group_split() %>%
  map(.f = function(datin) {
    panelplotted <- datin$Panel %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = GroupShort, y = Percent_NoStim,
          color = GroupShort, fill = GroupShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Phenograph Unstim by Group: ", panelplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/Phenograph_Unstim/",
      panelplotted, ".png"
    ),
    plotout,
    width = 8, height = 8
    )
  })

# stim by group
deltas_phenograph %>%
  group_by(Panel) %>%
  group_split() %>%
  map(.f = function(datin) {
    panelplotted <- datin$Panel %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = GroupShort, y = Percent_CMStim,
          color = GroupShort, fill = GroupShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Phenograph Stim by Group: ", panelplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/Phenograph_Stim/",
      panelplotted, ".png"
    ),
    plotout,
    width = 8, height = 8
    )
  })

# delta by group
deltas_phenograph %>%
  group_by(Panel) %>%
  group_split() %>%
  map(.f = function(datin) {
    panelplotted <- datin$Panel %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = GroupShort, y = Difference,
          color = GroupShort, fill = GroupShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Phenograph Delta by Group: ", panelplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/Phenograph_Delta/",
      panelplotted, ".png"
    ),
    plotout,
    width = 8, height = 8
    )
  })




# run beta regression models ===================================================

resultstable_unstim_lneg_lpos_phenograph <-
  run_unstim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                      correct_for_run = F, invert_OR_dir = F,
                      inputdeltadat = deltas_phenograph)
resultstable_stim_lneg_lpos_phenograph <-
  run_stim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                    correct_for_run = F, invert_OR_dir = F, correct_for_baseline = T,
                    inputdeltadat = deltas_phenograph)
resultstable_stim_lneg_lpos_phenograph_nbl <-
  run_stim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                    correct_for_run = F, invert_OR_dir = F, correct_for_baseline = F,
                    inputdeltadat = deltas_phenograph)




