# ==============================================================================
# 02_examine_distributions.R
# exploratory data analysis checking distributions of % positivity/lineage freq.
# ==============================================================================




# distribution of % outcomes ===================================================

# (arbitrary markers across varied mean values for examination)
# check unstim (can replace "Percent_NoStim" with "Percent_CMStim" for stim)
ggplot(
  deltas_final %>%
    filter(Panel == "APC" & CellType == "Mono" &
      Marker %in% c("IL1β", "IL10", "IL8", "PDL1")),
  aes(Difference, Marker,
      color = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.4, stat = "binline",
                      bins = 20,draw_baseline = F) +
  geom_vline(xintercept = 0, lty = 3) +
  theme_classic() +
  scale_x_continuous("Unstimulated %",
                     breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_y_discrete("Marker (Total Monocytes)", expand = c(0, 0)) +
  ggtitle("Example Data: Unstim") +
  scale_color_manual(values = palette_cohorts) +
  scale_fill_manual(values = palette_cohorts)


# deltas := stim - unstim
ggplot(
  deltas_final %>%
    filter(Panel == "APC" & CellType == "Mono" &
      Marker %in% c("IL1β", "IL10", "IL8", "PDL1")),
  aes(Difference, Marker,
    color = Group, fill = Group
  )
) +
  geom_density_ridges(
    alpha = 0.4,
    jittered_points = T, scale = 0.95, rel_min_height = 0.02,
    point_shape = "|", point_size = 2, size = 0.25,
    position = position_points_jitter(height = 0)
  ) +
  geom_vline(xintercept = 0, lty = 3) +
  theme_classic() +
  scale_x_continuous("Unstimulated %",
                     breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_y_discrete("Marker (Total Monocytes)", expand = c(0, 0)) +
  ggtitle("Example Data: Differences") +
  scale_color_manual(values = palette_cohorts) +
  scale_fill_manual(values = palette_cohorts)




# check corr btwn stim and % unstim ============================================
# should a true delta = stim - unstim be outcome (i.e., slope of one)?
# or should we instead model stim %s that are baseline adjusted?

plot_grid(
  deltas_final %>%
    filter(Panel == "APC" & CellType == "cDC1" & Marker == "IL10") %>%
    ggplot(data = ., aes(x = Percent_NoStim, y = Percent_CMStim)) +
    geom_point(alpha = 0.5, color = "darkred") +
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "lm", se = F, lty = 2, color = "darkred") +
    theme_few() +
    xlab("% No Stim") +
    ylab("% Stim") +
    ggtitle("APC // cDC1 // IL10"),
  deltas_final %>%
    filter(Panel == "APC" & CellType == "cDC2" & Marker == "IL27") %>%
    ggplot(data = ., aes(x = Percent_NoStim, y = Percent_CMStim)) +
    geom_point(alpha = 0.5, color = "darkred") +
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "lm", se = F, lty = 2, color = "darkred") +
    theme_few() +
    xlab("% No Stim") +
    ylab("% Stim") +
    ggtitle("APC // cDC2 // IL27")
)


# ideally looks like adjusting for baseline better than subtraction, even the
# latter subtraction delta method easier to understand
# test case -- beta regression fit
par(mfrow = c(2, 2))
deltas_final %>%
  filter(Panel == "APC" & CellType == "Mono" & Marker == "IL1β") %>%
  filter(Group != "BCG") %>%
  betareg(I(Percent_CMStim / 100) ~ Group + Percent_NoStim, data = .) %>%
  plot()




# plot features by run =========================================================
# checking for batch effects

c(
  "PlotsByRun_Unstim", "PlotsByRun_Stim", "PlotsByRun_Delta") %>%
  paste0("./FiguresTables/", current_date_label, "/", .) %>%
  sapply(dir.create)

# unstim by run
deltas_final %>%
  group_by(CellType) %>%
  group_split() %>%
  map(.f = function(datin) {
    celltypeplotted <- datin$CellType %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = RunShort, y = Percent_NoStim,
          color = RunShort, fill = RunShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Unstim by Run: ", celltypeplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/PlotsByRun_Unstim/",
      celltypeplotted, ".png"
    ),
    plotout,
    width = 5, height = 4.5
    )
  })

# stim by run
deltas_final %>%
  group_by(CellType) %>%
  group_split() %>%
  map(.f = function(datin) {
    celltypeplotted <- datin$CellType %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = RunShort, y = Percent_CMStim,
          color = RunShort, fill = RunShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Stim by Run: ", celltypeplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/PlotsByRun_Stim/",
      celltypeplotted, ".png"
    ),
    plotout,
    width = 5, height = 4.5
    )
  })

# delta by run
deltas_final %>%
  group_by(CellType) %>%
  group_split() %>%
  map(.f = function(datin) {
    celltypeplotted <- datin$CellType %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = RunShort, y = Difference,
          color = RunShort, fill = RunShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Delta by Run: ", celltypeplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/PlotsByRun_Delta/",
      celltypeplotted, ".png"
    ),
    plotout,
    width = 5, height = 4.5
    )
  })




# plot features by group =======================================================
# coarse, small plots to get fast impressions by celltype;
# in script 12, exports publication-ready figures

c(
  "", "PlotsByGroup_Unstim", "PlotsByGroup_Delta", "PlotsByGroup_Stim"
) %>%
  paste0("./FiguresTables/", current_date_label, "/", .) %>%
  sapply(dir.create)

# unstim by group
deltas_final %>%
  group_by(CellType) %>%
  group_split() %>%
  map(.f = function(datin) {
    celltypeplotted <- datin$CellType %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = GroupShort, y = Percent_NoStim,
          color = GroupShort, fill = GroupShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Unstim by Group: ", celltypeplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/PlotsByGroup_Unstim/",
      celltypeplotted, ".png"
    ),
    plotout,
    width = 5, height = 4.5
    )
  })

# stim by group
deltas_final %>%
  group_by(CellType) %>%
  group_split() %>%
  map(.f = function(datin) {
    celltypeplotted <- datin$CellType %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = GroupShort, y = Percent_CMStim,
          color = GroupShort, fill = GroupShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Stim by Group: ", celltypeplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/PlotsByGroup_Stim/",
      celltypeplotted, ".png"
    ),
    plotout,
    width = 5, height = 4.5
    )
  })

# delta by group
deltas_final %>%
  group_by(CellType) %>%
  group_split() %>%
  map(.f = function(datin) {
    celltypeplotted <- datin$CellType %>% first()
    plotout <-
      ggplot(
        datin,
        aes(
          x = GroupShort, y = Difference,
          color = GroupShort, fill = GroupShort
        )
      ) +
      geom_quasirandom(bandwidth = 0.25, size = 0.75, alpha = 0.75) +
      ggtitle(paste0("Delta by Group: ", celltypeplotted)) +
      facet_wrap(Marker ~ ., scales = "free_y") +
      theme_few() +
      theme(legend.position = "none")
    ggsave(paste0(
      "./FiguresTables/", current_date_label, "/PlotsByGroup_Delta/",
      celltypeplotted, ".png"
    ),
    plotout,
    width = 5, height = 4.5
    )
  })



# check associations with age ==================================================

# appear to be some immune associations with age,
# including some potential age*group effects: in minority of celltypes, but
# esp. in LTBI+/- vs BCG (probably not powered to examine/model these)
deltas_final %>%
  filter(Panel == "APC" & CellType == "cDC1" & Marker == "IL1β") %>%
  ggplot(data = ., aes(x = Age, y = Percent_NoStim, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, lty = 2) +
  theme_few() +
  scale_color_manual(values = palette_cohorts)

deltas_final %>%
  filter(Panel == "APC" & CellType == "Mono" & Marker == "PDL1") %>%
  ggplot(data = ., aes(x = Age, y = Percent_NoStim, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, lty = 2) +
  theme_few() +
  scale_color_manual(values = palette_cohorts)



