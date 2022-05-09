# ==============================================================================
# 09_correlation_heatmaps.R
# generate correlation heatmaps after covariate residualization
# ultimately excluded from final m.s., but useful for btwn marker correlation
# ==============================================================================




# fxn to format resids --> correlation =========================================

prep_residuals_for_heatmap <-
  function(resid_table, group_to_include) {
    resid_table <-
      resid_table %>%
      unite(col = "ID", c("Panel", "CellType", "Marker"), sep = " // ")

    PIDs_in_group <-
      deltas_final %>%
      filter(Group == group_to_include) %>%
      .$PID %>%
      unique()

    cor_table <-
      resid_table %>%
      dplyr::select(all_of(PIDs_in_group)) %>%
      as.matrix() %>%
      t() %>%
      set_colnames(resid_table$ID) %>%
      cor(method = "spearman", use = "complete")

    cor_table
  }

# BCG residuals could theoretically be from one of two models:
# LTBI- vs BCG or LTBI+ vs BCG models; just take average
harmonize_bcg_corrs <-
  function(tbl1, tbl2) {
    columns_to_include <-
      intersect(colnames(tbl1), colnames(tbl2))
    `+`(
      tbl1[columns_to_include, columns_to_include],
      tbl2[columns_to_include, columns_to_include]
    ) / 2
  }




# core plotting function =======================================================

dir.create(paste0("./FiguresTables/", current_date_label, "/CorHeatmap/"))

generate_heatmap <-
  function(resid_table = NULL, group_to_include = NULL,
           correlation_title, fileoutname, heatmap_input = NULL) {
    if (is.null(heatmap_input)) {

      # if no explicit residuals given as input, calculate
      heatmap_input <-
        prep_residuals_for_heatmap(resid_table, group_to_include)

      PIDs_in_group <-
        deltas_final %>%
        filter(Group == group_to_include) %>%
        .$PID %>%
        unique()

      # print sample size
      resid_table %>%
        dplyr::select(all_of(PIDs_in_group)) %>%
        as.matrix() %>%
        apply(., 2, function(x) {
          is.na(x) %>%
            any() %>%
            `!`()
        }) %>%
        sum() %>%
        paste0("N = ", .) %>%
        print()
    }

    # heatmap object
    heatmap_out <-
      Heatmap(
        column_title = correlation_title, column_title_gp = gpar(fontsize = 16),
        name = "Spearman Correlation",
        heatmap_input %>%
          set_colnames(., gsub("APC //|TCell //", "", colnames(.))) %>%
          set_rownames(., gsub("APC //|TCell//", "", colnames(.))),
        right_annotation =
          rowAnnotation(
            ` ` = colnames(heatmap_input) %>% grepl("APC", .) %>%
              if_else(., "APC", "TCell"),
            col = list(` ` = c(APC = "black", TCell = "grey")),
            simple_anno_size = unit(5, "point"),
            annotation_legend_param = list(
              title = "Panel",
              legend_direction = "horizontal",
              nrow = 1
            )
          ),
        bottom_annotation =
          columnAnnotation(
            ` ` = colnames(heatmap_input) %>% grepl("APC", .) %>%
              if_else(., "APC", "TCell"),
            col = list(` ` = c(APC = "black", TCell = "grey")),
            simple_anno_size = unit(5, "point"),
            show_legend = F
          ),
        row_names_gp = gpar(fontsize = 5.5),
        column_names_gp = gpar(fontsize = 5.5),
        col = colorRamp2(
          c(-1, -0.5, -0.15, 0, 0.15, 0.5, 1),
          c(
            "#053061", "#80b5cf", "#d3e3eb", "white",
            "#f2d1c2", "#f4a582", "#b2182b"
          ) %>% rev()
        ),
        heatmap_height = unit(7, "inches"), heatmap_width = unit(7, "inches"),
        show_column_dend = F,
        row_dend_side = "right",
        heatmap_legend_param = list(
          legend_height = unit(0.5, "cm"),
          legend_width = unit(6, "cm"),
          legend_direction = "horizontal"
        )
      )

    png(
      filename = paste0(
        "./FiguresTables/", current_date_label,
        "/CorHeatmap/", fileoutname, ".png"
      ),
      width = 8, height = 8, units = "in", res = 300
    )
    draw(heatmap_out,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom", merge_legend = T
    )
    dev.off()

    assign(x = fileoutname, value = heatmap_input, envir = .GlobalEnv)

    return(NULL)
  }




# unstimulated =================================================================

generate_heatmap(resultstable_unstim_lneg_lpos$Residuals,
  group_to_include = "LTBI-",
  correlation_title = "Correlations (unstimulated, LTBI-)",
  fileoutname = "unstim_ltbineg"
)

generate_heatmap(resultstable_unstim_lneg_lpos$Residuals,
  group_to_include = "LTBI+",
  correlation_title = "Correlations (unstimulated, LTBI+)",
  fileoutname = "unstim_ltbipos"
)

generate_heatmap(
  correlation_title = "Correlations (unstimulated, BCG)",
  fileoutname = "unstim_bcg",
  heatmap_input = harmonize_bcg_corrs(
    prep_residuals_for_heatmap(resultstable_unstim_lneg_bcg$Residuals, "BCG"),
    prep_residuals_for_heatmap(resultstable_unstim_lpos_bcg$Residuals, "BCG")
  )
)



# stimulated (bl-adjusted; memory) =============================================

generate_heatmap(resultstable_stim_lneg_lpos$Residuals,
  group_to_include = "LTBI-",
  correlation_title = "Correlations (Mtb-mem, LTBI-)",
  fileoutname = "stim_ltbineg"
)

generate_heatmap(resultstable_stim_lneg_lpos$Residuals,
  group_to_include = "LTBI+",
  correlation_title = "Correlations (Mtb-mem, LTBI+)",
  fileoutname = "stim_ltbipos"
)

generate_heatmap(
  correlation_title = "Correlations (Mtb-mem, BCG)",
  fileoutname = "stim_bcg",
  heatmap_input = harmonize_bcg_corrs(
    prep_residuals_for_heatmap(resultstable_stim_lneg_bcg$Residuals, "BCG"),
    prep_residuals_for_heatmap(resultstable_stim_lpos_bcg$Residuals, "BCG")
  )
)




# stimulated (nbl) =============================================================

generate_heatmap(
  resid_table = resultstable_stim_lneg_lpos_nbl$Residuals,
  group_to_include = "LTBI-",
  correlation_title = "Correlations (stimulated, nbl, LTBI-)",
  fileoutname = "stim_nbl_ltbineg"
)

generate_heatmap(resultstable_stim_lneg_lpos_nbl$Residuals,
  group_to_include = "LTBI+",
  correlation_title = "Correlations (stimulated, nbl, LTBI+)",
  fileoutname = "stim_nbl_ltbipos"
)

generate_heatmap(
  correlation_title = "Correlations (stimulated, nbl, BCG)",
  fileoutname = "stim_nbl_bcg",
  heatmap_input = harmonize_bcg_corrs(
    prep_residuals_for_heatmap(resultstable_stim_lneg_bcg_nbl$Residuals, "BCG"),
    prep_residuals_for_heatmap(resultstable_stim_lpos_bcg_nbl$Residuals, "BCG")
  )
)




