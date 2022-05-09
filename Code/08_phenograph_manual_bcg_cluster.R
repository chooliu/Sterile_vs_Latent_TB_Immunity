# ==============================================================================
# 08_phenograph_manual_bcg_cluster.R
# follow up on T-Cell TCluster26 in BCG participants
# ==============================================================================




# manual clustering of CD8+GMM+GranzB+ =========================================
# data from script 07 only compares LTBI- vs LTBI+
# additional manual clustering from EJ so we can examine in BCG participants

dat_deltas_phenograph_TCluster26 <-
  read_excel("./DataRaw/20220302/GMM CD8 GranzB for stats revised 3 3 22.xlsx") %>%
  mutate(Stim = if_else(grepl("CM", GROUPNAME), "CMStim", "NoStim")) %>%
  dplyr::select(PID, Stim, `CD8+GMM+GranzB+`) %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Marker", values_to = "Percent") %>%
  pivot_wider(., id_cols = c("PID", "Marker"), names_prefix = "Percent_",
              names_from = Stim, values_from = Percent) %>%
  mutate(Panel = "TCell", CellType = "Phenograph") %>%
  mutate(Difference = Percent_CMStim - Percent_NoStim,
         FreqOf = "Parent", Exclude_Marker = F, Exclude_CellType = F) %>%
  left_join(., metadata_consolated, by = c("PID") ) %>%
  filter(Group %in% c("BCG", "LTBI-", "LTBI+")) %>%
  mutate(Group = fct_drop(Group))


results_phenograph_TCluster26 <-
list(
  run_unstim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                      correct_for_run = F, invert_OR_dir = F,
                      inputdeltadat = dat_deltas_phenograph_TCluster26),
  run_unstim_analysis(c("LTBI-", "BCG"), "GroupLTBI-", "BCG", "LTBI-",
                      correct_for_run = F, invert_OR_dir = T,
                      inputdeltadat = dat_deltas_phenograph_TCluster26),
  run_unstim_analysis(c("LTBI+", "BCG"),
                      "GroupLTBI+", "BCG", "LTBI+",
                      correct_for_run = F, invert_OR_dir = T,
                      inputdeltadat = dat_deltas_phenograph_TCluster26)) %>%
  map_dfc(~ .x$Results %>% select(contains("P (LTBI"), contains("Odds ")))




# # stim analysis code would be similar (below) ================================
# however, this cluster doesn't pass default filtering criteria bc |deltas| low
# Mtb-stim doesn't change CD8+GMM+GranzB+ values substantially from unstim
run_stim_analysis(c("LTBI-", "LTBI+"),
                  "GroupLTBI+", "LTBI+", "LTBI-",
                  correct_for_run = F, invert_OR_dir = F, correct_for_baseline = T,
                  inputdeltadat = dat_deltas_phenograph_TCluster26)



