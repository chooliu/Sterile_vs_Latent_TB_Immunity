# ==============================================================================
# 03_multivariate_viz.R
# basic multivariate visualization: bar plots of parent freqs, PCA
# ==============================================================================




# bar plots ====================================================================

path_barplots <- paste0("./FiguresTables/", current_date_label, "/Barplots/")
dir.create(path_barplots)

ggplot(deltas_final %>%
         filter(FreqOf != "Parent" & Panel == "TCell" &
                  Exclude_Marker == F & Exclude_CellType == F) %>%
         arrange(Marker, -Percent_CMStim) %>%
         mutate(PID = as_factor(PID) %>% fct_inorder()),
       aes(PID, Percent_NoStim, fill = Marker)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_few(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90)) + #, legend.position = "none") +
  scale_y_continuous("% Unstim", expand = c(0, 0)) + 
  facet_grid( ~ Group, scales = "free", space = "free") +
  scale_fill_manual(values = qualitative_hcl(n = 12))
ggsave(paste0(path_barplots, "Parent_TCell.png"), width = 5, height = 4)

ggplot(deltas_final %>%
         filter(FreqOf != "Parent" & Panel == "APC" &
                  Exclude_Marker == F & Exclude_CellType == F) %>%
         arrange(Marker, -Percent_CMStim) %>%
         mutate(PID = as_factor(PID) %>% fct_inorder()),
       aes(PID, Percent_NoStim, fill = Marker)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_few(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90)) + #, legend.position = "none") +
  scale_y_continuous("% Unstim", expand = c(0, 0)) + 
  facet_grid( ~ Group, scales = "free", space = "free") +
  scale_fill_manual(values = qualitative_hcl(n = 5))
ggsave(paste0(path_barplots, "Parent_APC.png"), width = 5, height = 4)




# generic PCA fxn ==============================================================
# colored by group; can also +/- shape by run, PID labels

path_pcaplots <- paste0("./FiguresTables/", current_date_label, "/PCA/")
dir.create(path_pcaplots)

fxn_run_pca <-
  function(pca_input, title_in, file_out) {
    pca_soln <-
      pca_input %>%
      dplyr::select(-PID) %>%
      prcomp(scale = T)

    pca_df <-
      bind_cols(
        PID = pca_input$PID,
        PC = pca_soln %>%
          .$x %>%
          .[, 1:2] %>%
          as_tibble()
      ) %>%
      right_join(metadata_consolated, ., by = "PID")

    pca_varexp <-
      (pca_soln$sdev)^2 %>%
      `/`(., sum(.)) %>%
      `*`(100) %>%
      formatC(digits = 1, format = "f") %>%
      paste0("PC", 1:length(.), " (", ., "%)")

    plotout <-
      ggplot(
        pca_df,
        aes(PC1, PC2, fill = Group, color = Group)
      ) +
      geom_point(aes(shape = Run), size = 2) +
      theme_few() +
      geom_text_repel(aes(label = PID)) +
      geom_convexhull(alpha = 0.2) +
      scale_color_manual(values = palette_cohorts) +
      scale_fill_manual(values = palette_cohorts) +
      scale_shape_manual(values = palette_run_shapes) +
      xlab(pca_varexp[1]) +
      xlab(pca_varexp[2]) +
      ggtitle(title_in)

    ggsave(paste0(path_pcaplots, file_out, ".png"), plotout,
      width = 5, height = 4
    )
    plotout
  }




# unstim PCAs ==================================================================

# parent populations -----------------------------------------------------------
# id "10146-01" & "10378-01" potential APC outliers,
# but no external rationale for exclusion

pca_input_parent_unstim <-
  deltas_final %>%
  filter(FreqOf != "Parent" &
    Exclude_Marker == F & Exclude_CellType == F) %>%
  mutate(Outcome = Percent_NoStim / 100) %>%
  mutate(Outcome = if_else(Outcome == 1, Outcome - 1e-6, Outcome)) %>%
  mutate(Outcome = if_else(Outcome == 0, Outcome + 1e-6, Outcome)) %>%
  mutate(OutcomeForPCA = log((Outcome + 1e-6) / (1 - Outcome + 1e-6))) %>%
  pivot_wider(
    id_cols = PID,
    names_from = c("CellType", "Marker"),
    values_from = OutcomeForPCA
  )

fxn_run_pca(
  pca_input_parent_unstim %>% .[, grepl("PID|APC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ],
  "APC Parent Frequencies, % Unstim", "APC_Parent"
)
fxn_run_pca(
  pca_input_parent_unstim %>% .[, grepl("PID|PBMC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ],
  "PBMC Parent Frequencies, % Unstim", "TCell_Parent"
)

# percent positivities ---------------------------------------------------------

pca_input_markers_unstim <-
  deltas_final %>%
  filter(FreqOf == "Parent" &
    Exclude_Marker == F & Exclude_CellType == F) %>%
  mutate(Outcome = Percent_NoStim / 100) %>%
  mutate(Outcome = if_else(Outcome == 1, Outcome - 1e-6, Outcome)) %>%
  mutate(Outcome = if_else(Outcome == 0, Outcome + 1e-6, Outcome)) %>%
  mutate(OutcomeForPCA = log((Outcome) / (1 - Outcome))) %>%
  pivot_wider(
    id_cols = PID,
    names_from = c("Panel", "CellType", "Marker"),
    values_from = OutcomeForPCA
  )

fxn_run_pca(
  pca_input_markers_unstim %>% .[, grepl("PID|APC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ] %>%
    .[, apply(., 2, function(x) {
      sd(x, na.rm = T) %>%
        `==`(0, .) %>%
        `!`()
    }) %>%
      (function(x) {
        is.na(x) | x == T
      })],
  "APC Markers, % Unstim", "APC_Markers"
)
fxn_run_pca(
  pca_input_markers_unstim %>% .[, grepl("PID|TCell", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ] %>%
    .[, apply(., 2, function(x) {
      sd(x, na.rm = T) %>%
        `==`(0, .) %>%
        `!`()
    }) %>%
      (function(x) {
        is.na(x) | x == T
      })],
  "PBMC Markers, % Unstim", "TCell_Markers"
)




# stim PCAs ====================================================================

# parent populations -----------------------------------------------------------

pca_input_parent_stim <-
  deltas_final %>%
  filter(FreqOf != "Parent" &
    Exclude_Marker == F & Exclude_CellType == F) %>%
  mutate(Outcome = Percent_CMStim / 100) %>%
  mutate(Outcome = if_else(Outcome == 1, Outcome - 1e-6, Outcome)) %>%
  mutate(Outcome = if_else(Outcome == 0, Outcome + 1e-6, Outcome)) %>%
  mutate(OutcomeForPCA = log((Outcome + 1e-6) / (1 - Outcome + 1e-6))) %>%
  pivot_wider(
    id_cols = PID,
    names_from = c("CellType", "Marker"),
    values_from = OutcomeForPCA
  ) %>%
  .[apply(., 1, function(x) {
    is.na(x) %>%
      any() %>%
      `!`()
  }), ]

fxn_run_pca(
  pca_input_parent_stim %>% .[, grepl("PID|APC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ],
  "APC Parent Frequencies, % Stim", "APC_Parent_Stim"
)
fxn_run_pca(
  pca_input_parent_stim %>% .[, grepl("PID|PBMC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ],
  "PBMC Parent Frequencies, % Stim", "TCell_Parent_Stim"
)

# percent positivities ---------------------------------------------------------

pca_input_markers_stim <-
  deltas_final %>%
  filter(FreqOf == "Parent" &
    Exclude_Marker == F & Exclude_CellType == F) %>%
  mutate(Outcome = Percent_CMStim / 100) %>%
  mutate(Outcome = if_else(Outcome == 1, Outcome - 1e-6, Outcome)) %>%
  mutate(Outcome = if_else(Outcome == 0, Outcome + 1e-6, Outcome)) %>%
  mutate(OutcomeForPCA = log((Outcome + 1e-6) / (1 - Outcome + 1e-6))) %>%
  pivot_wider(
    id_cols = PID,
    names_from = c("Panel", "CellType", "Marker"),
    values_from = OutcomeForPCA
  )

fxn_run_pca(
  pca_input_markers_stim %>% .[, grepl("PID|APC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ] %>%
    .[, apply(., 2, function(x) {
      sd(x, na.rm = T) %>%
        `==`(0, .) %>%
        `!`()
    }) %>%
      (function(x) {
        is.na(x) | x == T
      })],
  "APC Markers, % Stim", "APC_Markers_Stim"
)
fxn_run_pca(
  pca_input_markers_stim %>% .[, grepl("PID|TCell", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ] %>%
    .[, apply(., 2, function(x) {
      sd(x, na.rm = T) %>%
        `==`(0, .) %>%
        `!`()
    }) %>%
      (function(x) {
        is.na(x) | x == T
      })],
  "PBMC Markers, % Stim", "TCell_Markers_Stim"
)




# delta PCAs ===================================================================

# parent populations -----------------------------------------------------------

pca_input_parent_delta <-
  deltas_final %>%
  filter(FreqOf != "Parent" &
    Exclude_Marker == F & Exclude_CellType == F) %>%
  mutate(OutcomeForPCA = Percent_CMStim - Percent_NoStim) %>%
  pivot_wider(
    id_cols = PID,
    names_from = c("CellType", "Marker"),
    values_from = OutcomeForPCA
  ) %>%
  .[apply(., 1, function(x) {
    is.na(x) %>%
      any() %>%
      `!`()
  }), ]

fxn_run_pca(
  pca_input_parent_delta %>% .[, grepl("PID|APC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ],
  "APC Parent Frequencies, Deltas", "APC_Parent_Delta"
)
fxn_run_pca(
  pca_input_parent_delta %>% .[, grepl("PID|PBMC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ],
  "PBMC Parent Frequencies, Deltas", "TCell_Parent_Delta"
)

# percent positivities ---------------------------------------------------------

pca_input_markers_delta <-
  deltas_final %>%
  filter(FreqOf == "Parent" &
    Exclude_Marker == F & Exclude_CellType == F) %>%
  mutate(OutcomeForPCA = Percent_CMStim - Percent_NoStim) %>%
  pivot_wider(
    id_cols = PID,
    names_from = c("Panel", "CellType", "Marker"),
    values_from = OutcomeForPCA
  )

fxn_run_pca(
  pca_input_markers_delta %>% .[, grepl("PID|APC", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ] %>%
    .[, apply(., 2, function(x) {
      sd(x, na.rm = T) %>%
        `==`(0, .) %>%
        `!`()
    }) %>%
      (function(x) {
        is.na(x) | x == T
      })],
  "APC Markers, Stim Deltas", "APC_Markers_Delta"
)
fxn_run_pca(
  pca_input_markers_delta %>% .[, grepl("PID|TCell", names(.))] %>%
    .[apply(., 1, function(x) {
      is.na(x) %>%
        any() %>%
        `!`()
    }), ] %>%
    .[, apply(., 2, function(x) {
      sd(x, na.rm = T) %>%
        `==`(0, .) %>%
        `!`()
    }) %>%
      (function(x) {
        is.na(x) | x == T
      })],
  "PBMC Markers, Stim Deltas", "TCell_Markers_Delta"
)



