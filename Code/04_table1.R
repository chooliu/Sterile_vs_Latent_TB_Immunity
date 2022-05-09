# ==============================================================================
# 04_table1.R
# participant demographic features by cohort
# ==============================================================================




# check distribution of metadata, by group -------------------------------------
ggplot(metadata_consolated %>%
         filter(Group %in% c("BCG", "LTBI-", "LTBI+") &
                  Run != "Run4"),
       aes(x = Age,  # examine also BMI, ExposureScore
           y = Group, color = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.4,
                      jittered_points = T, scale = 0.95, rel_min_height = 0.02,
                      point_shape = "|", point_size = 2, size = 0.25,
                      position = position_points_jitter(height = 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(values = palette_cohorts) +
  scale_fill_manual(values = palette_cohorts)  +
  theme_classic() + theme(legend.position = "none")



# core table 1 variable summaries ----------------------------------------------

table1_allsubj <-
  deltas_final %>%
  ungroup %>%
  filter(!duplicated(PID)) %>%
  bind_rows(., mutate(., Group = "All")) %>%
  group_by(Group) %>%
  summarize(N = n(),
            Age = summaryMedianRange(Age),
            `Sex (Female)` = summaryCountPercent(Sex, "Female"),
            `BMI` = summaryMedianRange(BMI, na.rm = T),
            `Exposure Score` = summaryMedianRange(ExposureScore, na.rm = T),
            `BCG Scar (present)` = summaryCountPercent(BCGScar, "Present"),
            `Race (South Asian)` = summaryCountPercent(Race, "South Asian"),
            `Race (Asian)` = summaryCountPercent(Race, "A"),
            `Race (Black)` = summaryCountPercent(Race, "B"),
            `Race (Hispanic White)` = summaryCountPercent(Race, "H/W"),
            `Race (Non-Hispanic White)` = summaryCountPercent(Race, "W"),
            `Run 1` = summaryCountPercent(Run, "Run1"),
            `Run 2` = summaryCountPercent(Run, "Run2"),
            `Run 3` = summaryCountPercent(Run, "Run3"),
            `Run 4` = summaryCountPercent(Run, "Run4") # later omitted due to batch effects
            ) %>%
  t() %>%
  .[ , c(2:4, 1)]

# p-values ---------------------------------------------------------------------
# note: previously separated by panel (APC, T-Cell) since not all particip.
# have both, final paper presents the the all subjects version

input_table1_ptests <-
  deltas_final %>%
  mutate(Group = fct_drop(Group)) %>%
  bind_rows(., mutate(., Panel = "All")) %>%
  group_by(Panel) %>%
  group_split() %>%
  map( ~ filter(.x, !duplicated(PID)))

input_table1_ptests_nobcg <-
  map(input_table1_ptests,
      ~ filter(.x, Group != "BCG")) %>%
  map(~ mutate(.x, Group = fct_drop(Group)))


get_tbl1_pvals <-
  function(data_by_panel) {
  list(
    
    age_pval =
         data_by_panel %>%
         map(~ anova(lm(Age ~ Group, .x),
                     lm(Age ~ 1, .x)) ) %>%
       map_dbl( ~ .x$`Pr(>F)`[2]),
    
     sex_pval =
      data_by_panel %>%
       map(~ filter(.x, !is.na(.x$Sex))) %>%
       map(~ mutate(.x, Group = fct_drop(Group))) %>%
       map(~ table(.x$Sex, .x$Group) %>% fisher.test()) %>%
      map_dbl( ~ .x$p.value),

    bmi_pval =
      data_by_panel %>%
      map(~ anova(lm(BMI ~ Group, .x),
                  lm(BMI ~ 1, .x)) ) %>%
      map_dbl( ~ .x$`Pr(>F)`[2]),

    
    exposure_pval =
      data_by_panel %>%
      map(~ anova(lm(ExposureScore ~ Group, .x),
                  lm(ExposureScore ~ 1, .x)) ) %>%
      map_dbl( ~ .x$`Pr(>F)`[2]),
    
    bcgscar_pval =
      data_by_panel %>%
      map(~ filter(.x, !is.na(.x$BCGScar))) %>%
      map(~ mutate(.x, Group = fct_drop(Group))) %>%
      map(~ table(.x$BCGScar, .x$Group) %>% fisher.test()) %>%
      map_dbl( ~ .x$p.value),
    
    race_pval =
        data_by_panel %>%
          map(~ filter(.x, !is.na(.x$Race))) %>%
          map(~ mutate(.x, Group = fct_drop(Group))) %>%
          map_dbl(function(x) {
            if ("BCG" %in% x$Group) { # error if only LTBI+/- (exclusively S.A.)
              table(x$Race, x$Group) %>% fisher.test() %>% .$p.value
              } else {
                NA_real_
              }}) 
          )
    }


# final table ------------------------------------------------------------------
# (replace with ~.x[1] or [2] to extract p-vals for APC or T-Cell subjects only)

pvalues_tbl1_threegroup <- get_tbl1_pvals(input_table1_ptests)
pvalues_tbl1_twogroup <- get_tbl1_pvals(input_table1_ptests_nobcg)

table1_allsubj <-
  cbind(table1_allsubj,
        `P (three group)` =
          c(NA_integer_, NA_integer_, map_dbl(pvalues_tbl1_threegroup, ~ .x[3]),
            rep(NA_integer_, 8)) %>% formatC(digits = 2),
        `P (LTBI- vs +)` =
          c(NA_integer_, NA_integer_, map_dbl(pvalues_tbl1_twogroup, ~ .x[3]),
            rep(NA_integer_, 8)) %>% formatC(digits = 2))



