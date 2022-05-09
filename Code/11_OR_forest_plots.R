# ==============================================================================
# 11_OR_forest_plots.R
# generate forest-like plots for OR (95% CI)
# ==============================================================================




# core plotting functions ======================================================
# note: each is marginally modified in a .svg editor afterwards
# (making spacing between OR and FDR more compact, adding arrows  + text header)

# error-bars / 'trees' to show confidence intervals
forest_plot_OR_panel <-
  function(resultstable,
           xaxis_title = "", yaxis_label,
           xaxis_labels,
           txt_positions,
           filter_signif = T,
           sort_style_marker = T,
           plot_color_scheme = palette_cohorts) {
    
    position_or <- txt_positions[1]
    position_grp1 <- txt_positions[2]
    position_grp2 <- txt_positions[3]
    position_fdr <- txt_positions[4]
    
    # extract key columns
    # (previously showed predicted means, too cluttered for final figure)
    plot_input <-
      resultstable %>%
      select(Panel, CellType, Marker, contains("Signif"),
             contains("Raw"), contains("Odds R"), contains("FDR")) %>%
      set_names(c("Panel", "CellType", "Marker", "Signif",
                  "Pred1", "Pred2", "OR", "FDR"))
    
    # show signif only?
    filter_markers <- T
    if (filter_signif) {
      filter_markers <- plot_input$Signif == "*"
    }
    
    # sort labels in style for plotting activation markers
    # (celltype --> protein name)
    if (sort_style_marker) {
      plot_input <-
        plot_input[filter_markers, ] %>%
        mutate(CellType = str_pad(CellType, width = 16, "right", pad = " ")) %>%
        mutate(DCs = !grepl("cDC|pDC|dnDC|pDC", Marker, ignore.case = F)) %>%
        arrange(Panel, DCs, CellType, Marker) %>%
        mutate(CellType = factor(CellType) %>% fct_inorder,
               Marker = factor(Marker))
        plot_input$Marker <- fct_rev(plot_input$Marker)
      } else {
        plot_input <-
          plot_input[filter_markers, ] %>%
          mutate(CellType = str_pad(CellType, width = 16, "right", pad = " ")) %>%
          mutate(DCs = !grepl("cDC|pDC|dnDC|pDC", Marker, ignore.case = F)) %>%
          arrange(Panel, DCs, Marker) %>%
          mutate(CellType = factor(CellType),
                 Marker = factor(Marker) %>% fct_inorder)
        plot_input$Marker <- fct_rev(plot_input$Marker)
      }
    
    # formatting text / extracting ORs
    plot_input <-
      plot_input %>%
      mutate(FDR = if_else(grepl("e-", FDR), gsub("e-", "%*%10^{-", FDR) %>%
                             paste0("", ., "}"), FDR) %>%
               paste0(., if_else(Signif == "*", "", "~(ns)")) %>%
               as.expression %>% as.character,
             OR = OR %>% gsub("\n ", " ", .),
             OR_color = if_else(OR <= 1, "LTBI-", "LTBI+"),
             OR_edit = OR %>% gsub("\n\\(| - |\\)", "_", .)) %>%
        separate(col = OR_edit, into = c("mid", "low", "high"),
                 sep = "_", convert = T, remove = F) %>%
      mutate_at(.vars = c("mid", "low", "high"),
                .funs = log)
    
    # alternate grey/white by celltype
      plot_background <-
        plot_input %>%
        select(Panel, CellType) %>%
        mutate(plot_rect_color = CellType %in%
                 (plot_input$CellType %>% unique %>% .[c(T, F)])) %>%
        filter(!duplicated(CellType))
      
      # final plot
      ggplotout <-
        ggplot(data = plot_input) + 
        geom_rect(data = plot_background,
                  aes(fill = plot_rect_color),
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 1) +
        geom_errorbar(aes(x = Marker, ymin = low, ymax = high,
                          color = OR_color), width = 0.7) +
        geom_hline(yintercept = c(min(xaxis_labels),
                                  max(xaxis_labels)), lty = 1) +
        geom_hline(yintercept = max(xaxis_labels) + 0.1,
                   lty = 1, color = "#ffffff", size = 2) +
        geom_hline(yintercept = 0, lty = 3) +
        geom_point(aes(x = Marker, y = mid, color = OR_color),
                   shape = 19, size = 1.25) +
        scale_x_discrete(yaxis_label) + 
        scale_y_continuous(
          xaxis_title,
          breaks = xaxis_labels, labels = xaxis_labels, 
          limits = c(min(xaxis_labels), max(xaxis_labels) + 5),
          expand = c(0, 0)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            axis.text.x = element_text(hjust = 0, vjust = 0, size = 10),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
            axis.title = element_text(hjust = 0, size = 11),
            axis.title.y = element_text(angle = 0, hjust = 0, vjust = 1, size = 13),
            strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1, margin = margin(b = 3)),
            strip.placement = "outside", panel.spacing = unit(0.075, "cm"),
            axis.ticks.x = element_line(color = "black", size = 0.5)) +
        scale_fill_manual(values = c("white", "#F0F0F0")) +
        scale_alpha_manual(values = c(1, 0)) +
        scale_color_manual(values = plot_color_scheme) + 
        facet_grid(CellType ~ ., scale = "free_y", space = "free_y", switch = "y") +
        coord_flip() +
        geom_text(aes(x = Marker, y = position_or,
                      label = formatC(exp(mid), format = "f", digits = 2)),
                  hjust = 0, size = 3.5) +
        geom_text(aes(x = Marker, y = position_grp1,
                      label = formatC(as.numeric(Pred1), format = "f", digits = 1)),
                  hjust = 0, size = 3.5) +
        geom_text(aes(x = Marker, y = position_grp2,
                      label = formatC(as.numeric(Pred2), format = "f", digits = 1)),
                  hjust = 0, size = 3.5) +
        geom_text(aes(x = Marker, y = position_fdr, label = FDR),
                  hjust = 0, size = 3.5, parse = T)
      
      ggplotout
  }

# export folder
dir.create(paste0("./FiguresTables/", current_date_label, "/Forest/"))

# for labeling on left hand side in big 3-panel figures
hypothesis_labels <-
  c("Unstimulated", "Mtb-Stimulated", "Mtb-Memory") %>% 
  str_pad(., width = 20, side = "right")




# LTBI- vs LTBI+ ===============================================================

xlabels_LvL <- c(-3, -2, -1, 0, 1, 2)
txtpositions_LvL <- c(2.5, 4.5, 3.5, 5.5) 

# unstim
forest_LvL_unstim <-
    forest_plot_OR_panel(resultstable_unstim_lneg_lpos$Results,
                         xaxis_title = "ln(Odds Ratio): LTBI/Resister",
                         xaxis_labels = xlabels_LvL,
                         txt_positions = txtpositions_LvL,
                         yaxis_label = "")
ggsave(forest_LvL_unstim, width = 6, height = 3,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/LvL_unstim.pdf"),
       device = cairo_pdf)
paste0("open ./FiguresTables/20220325/Forest/LvL_unstim.pdf") %>% system

# Mtb-stim (no-baseline adj)
forest_LvL_stim_nbl <-
  forest_plot_OR_panel(resultstable_stim_lneg_lpos_nbl$Results %>% select(!contains("Delta")),
                       xaxis_title = "ln(Odds Ratio): LTBI/Resister",
                       xaxis_labels = xlabels_LvL,
                       txt_positions = txtpositions_LvL,
                       yaxis_label = "")

ggsave(forest_LvL_stim_nbl, width = 6, height = 5,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/LvL_stim_nbl.pdf"),
       device = cairo_pdf)
paste0("open ./FiguresTables/20220325/Forest/LvL_stim_nbl.pdf") %>% system

# Mtb-memory
forest_LvL_delta <-
  forest_plot_OR_panel(resultstable_stim_lneg_lpos$Results %>% select(!contains("Stim %")),
                     xaxis_title = "ln(Odds Ratio): LTBI/Resister",
                     xaxis_labels = xlabels_LvL,
                     txt_positions = txtpositions_LvL,
                     yaxis_label = "")
ggsave(forest_LvL_delta, width = 6, height = 2,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/LvL_delta.pdf"),
       device = cairo_pdf)
paste0("open ./FiguresTables/20220325/Forest/LvL_delta.pdf") %>% system

# parent populations
xlabels_LvLparent <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
txtpositions_LvLparent <- c(2, 3.75, 3, 4.5)
forest_LvLparent <-
  forest_plot_OR_panel(resultstable_unstim_lneg_lpos_parent$Results,
                       xaxis_title = "ln(Odds Ratio): LTBI/Resister",
                     xaxis_labels = xlabels_LvLparent,
                     txt_positions = txtpositions_LvLparent,
                     yaxis_label = "",
                     filter_signif = F)
ggsave(forest_LvLparent, width = 7, height = 4,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/LvL_parent.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/LvL_parent.pdf'") %>% system

# phenograph
xlabels_LvLpheno <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
txtpositions_LvLpheno <- c(2, 3.5, 2.75, 4.25)
forest_LvLpheno <-
  plot_grid(
    forest_plot_OR_panel(resultstable_unstim_lneg_lpos_phenograph$Results %>%
                           mutate(CellType = Panel),
                         xaxis_title = "",
                         xaxis_labels = xlabels_LvLpheno,
                         txt_positions = txtpositions_LvLpheno,
                         yaxis_label = hypothesis_labels[1]),
    forest_plot_OR_panel(resultstable_stim_lneg_lpos_phenograph_nbl$Results %>%
                           mutate(CellType = Panel) %>% select(!contains("Delta")),
                         xaxis_title = "",
                         xaxis_labels = xlabels_LvLpheno,
                         txt_positions = txtpositions_LvLpheno,
                         yaxis_label = hypothesis_labels[2]),
    forest_plot_OR_panel(resultstable_stim_lneg_lpos_phenograph$Results %>%
                           mutate(CellType = Panel) %>% select(!contains("Stim %")),
                         xaxis_title = "ln(Odds Ratio): LTBI/Resister",
                         xaxis_labels = xlabels_LvLpheno,
                         txt_positions = txtpositions_LvLpheno,
                         yaxis_label = hypothesis_labels[3]),
    nrow = 3, rel_heights = c(3, 1.75, 3), align = "h")
ggsave(forest_LvLpheno, width = 8, height = 6,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/LvL_pheno.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/LvL_pheno.pdf'") %>% system




# BCG vs LTBI- =================================================================

xlabels_BvLneg <- seq(-2.5, 2, 0.5)
txtpositions_BvLneg <- c(2.5, 4.5, 3.5, 5.5)

# unstim
fig_BvLneg_unstim <-
  forest_plot_OR_panel(resultstable_unstim_lneg_bcg$Results,
                       xaxis_title = "ln(Odds Ratio): Resister/BCG",
                       xaxis_labels = xlabels_BvLneg,
                       txt_positions = txtpositions_BvLneg,
                       yaxis_label = "",
                       plot_color_scheme = palette_cohorts[c(1, 2)] %>% set_names(NULL))
ggsave(fig_BvLneg_unstim, width = 8, height = 8,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/BvLneg_unstim.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/BvLneg_unstim.pdf'") %>% system

# Mtb-stim
fig_BvLneg_stim_nbl <-
  forest_plot_OR_panel(resultstable_stim_lneg_bcg_nbl$Results %>% select(!contains("Delta %")),
                       xaxis_title = "ln(Odds Ratio): Resister/BCG",
                       xaxis_labels = xlabels_BvLneg,
                       txt_positions = txtpositions_BvLneg,
                       yaxis_label = "",
                       plot_color_scheme = palette_cohorts[c(1, 2)] %>% set_names(NULL))
ggsave(fig_BvLneg_stim_nbl, width = 8, height = 8,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/BvLneg_stim_nbl.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/BvLneg_stim_nbl.pdf'") %>% system

# Mtb-memory
forest_BvLneg_delta <-
  forest_plot_OR_panel(resultstable_stim_lneg_bcg$Results %>% select(!contains("Stim %")),
                       xaxis_title = "ln(Odds Ratio): Resister/BCG",
                       xaxis_labels = xlabels_BvLneg,
                       txt_positions = txtpositions_BvLneg,
                       yaxis_label = "",
                       plot_color_scheme = palette_cohorts[c(1, 2)] %>% set_names(NULL))
ggsave(forest_BvLneg_delta, width = 8, height = 8,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/BvLneg_delta.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/BvLneg_delta.pdf'") %>% system




# BCG vs LTBI+ =================================================================

xlabels_BvLpos <- seq(-3, 4, 1)
txtpositions_BvLpos <- c(2.5, 4.5, 3.5, 5.5) + 2

# unstim
fig_BvLpos_unstim <-
  forest_plot_OR_panel(resultstable_unstim_lpos_bcg$Results,
                       xaxis_title = "ln(Odds Ratio): LTBI/BCG",
                       xaxis_labels = xlabels_BvLpos,
                       txt_positions = txtpositions_BvLpos,
                       yaxis_label = "",
                       plot_color_scheme = palette_cohorts[c(1, 3)] %>% set_names(NULL))
ggsave(fig_BvLpos_unstim, width = 8, height = 9,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/BvLpos_unstim.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/BvLpos_unstim.pdf'") %>% system

# Mtb-stim
fig_BvLpos_stim_nbl <-
  forest_plot_OR_panel(resultstable_stim_lpos_bcg_nbl$Results %>% select(!contains("Delta %")),
                       xaxis_title = "ln(Odds Ratio): LTBI/BCG",
                       xaxis_labels = xlabels_BvLpos,
                       txt_positions = txtpositions_BvLpos,
                       yaxis_label = "",
                       plot_color_scheme = palette_cohorts[c(1, 3)] %>% set_names(NULL))
ggsave(fig_BvLpos_stim_nbl, width = 8, height = 6,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/BvLpos_stim_nbl.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/BvLpos_stim_nbl.pdf'") %>% system

# Mtb-memory
forest_BvLpos_delta <-
  forest_plot_OR_panel(resultstable_stim_lpos_bcg$Results %>% select(!contains("Stim %")),
                       xaxis_title = "ln(Odds Ratio): LTBI/BCG",
                       xaxis_labels = xlabels_BvLpos,
                       txt_positions = txtpositions_BvLpos,
                       yaxis_label = "",
                       plot_color_scheme = palette_cohorts[c(1, 3)] %>% set_names(NULL))
ggsave(forest_BvLpos_delta, width = 8, height = 5.5,
       filename = paste0("./FiguresTables/", current_date_label, "/Forest/BvLpos_delta.pdf"),
       device = cairo_pdf)
paste0("open './FiguresTables/20220325/Forest/BvLpos_delta.pdf'") %>% system




