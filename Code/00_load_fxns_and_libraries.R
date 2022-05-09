# ==============================================================================
# 00_load_fxns_and_libraries.R
# basic setup, R packages & custom functions
# ==============================================================================




# libraries ====================================================================

# data manipulation
library(magrittr)
library(tidyverse)
library(readxl)

# statistical testing
library(betareg)
library(lmtest)

# visualization & presentation
library(ComplexHeatmap)
library(cowplot)
library(circlize)
library(ggridges)
library(ggbeeswarm)
library(ggthemes)
library(writexl)
library(colorspace)
library(ggrepel)
library(ggConvexHull)
library(scales)

# export for version control
writeLines(
  capture.output(sessionInfo()),
  paste0("package_versions_", current_date_label, ".txt"))




# palettes for plotting ========================================================

palette_cohorts <-
  c(BCG = "#762a83", `LTBI-` = "#1f78b4", `LTBI+` = "#e41a1c")
palette_cohortsshort <-
  c(B = "#762a83", `L-` = "#1f78b4", `L+` = "#e41a1c")

palette_cohorts_shapes <-
  c(BCG = 15, `LTBI-` = 1, `LTBI+` = 19)
palette_run_shapes <-
  c(Run1 = 1, Run3 = 4, Run2 = 19, Run4 = 2)




# custom functions =============================================================
# for loading percent positivity FlowJo output
# (somewhat complex function i've drafted for Weinberg lab projects,
# typicallytransferred via a standard wide format .xlsx file)



# extract_flowcy_from_single_xlsx()
# ------------------------------------------------------------------------------
# input arguments:
# - filepath: location of .xlsx file of interest
# - marker_name_str: string to remove from markers (e.g., "Freq. of Parent")
#       (optional, could just replace with "" for no removal)
#       currently assumes that marker columsn have term "Singlets" in them
# - sample metadata tab: if patient ID / sample info tab exists in .xlsx
#       (if user doesn't supply, R will print a list of tabs and ask for input)
# output:
# - a list with three items, accessible using
#   * list_name$Percents: marker percent positivities in long (tidy) format
#   * list_name$Counts: raw event counts in long (tidy) format
#   * list_name$Metadata: metadata tab, read in its entirety
# ------------------------------------------------------------------------------

extract_flowcy_from_single_xlsx <-
  function(filepath,
           marker_name_str = c("Singlets/Live/Size/Cleanup/dump-/| \\| Freq. of Parent| \\| Freq. of DR+| \\ Count| \\| Freq. of Cleanup|Time/Singlets/Live/Size/Cleanup/"),
           metadata_tab_num = NULL) {
    
    paste0("reading in '", filepath, "'\n") %>% cat()
    tab_names <- excel_sheets(filepath)
    tab_names %>% print
    
    if (is.null(metadata_tab_num)) {
      paste0(seq(length(tab_names)), " ... ", tab_names, "\n") %>% cat()
      metadata_tab_num <-
        readline(prompt = "Are any of these tabs sample metadata? (input # or 'no'):  ") %>%
        as.numeric()
    }

    sheet_nums_to_read <-
      if (is.na(metadata_tab_num)) {
        1:length(tab_names)
      } else {
        setdiff(1:length(tab_names), metadata_tab_num)
      }

    percents <-
      map(
        sheet_nums_to_read,
        function(tab) {
          
          excel_file <- read_excel(filepath, sheet = tab_names[tab])
          columns_with_percents <- names(excel_file) %>% .[grepl("Singlets", .) & grepl("Freq.", .)]
          columns_with_counts <- names(excel_file) %>% .[grepl("Singlets", .) & grepl("Count", .)]
          columns_metadata <- names(excel_file) %>%
            setdiff(columns_with_percents) %>% setdiff(columns_with_counts) 
          
          if (length(columns_with_percents) == 0) { return(NULL) }
          
          cleaned_dat_PERCENT <-
            excel_file %>%
            select(-all_of(columns_with_counts)) %>%
            pivot_longer(
              cols = all_of(columns_with_percents),
              names_to = "Marker_Full", values_to = "Percent"
            ) %>%
            mutate(
              CellType = tab_names[tab],
              SourceFile = filepath,
            ) %>%
            select(all_of(c(setdiff(names(.), "Percent"), "Percent"))) # *
          # * convoluted, but makes percent last column for tidiness
        }
      ) %>%
      bind_rows()

    
    counts <-
      map(
        sheet_nums_to_read,
        function(tab) {
          
          excel_file <- read_excel(filepath, sheet = tab_names[tab])
          columns_with_percents <- names(excel_file) %>% .[grepl("Singlets", .) & grepl("Freq.", .)]
          columns_with_counts <- names(excel_file) %>% .[grepl("Singlets", .) & grepl("Count", .)]
          columns_metadata <-
            names(excel_file) %>%
            setdiff(columns_with_percents) %>% setdiff(columns_with_counts) 
          
          if (length(columns_with_counts) == 0) { return(NULL) }
          
          cleaned_dat_COUNTS <-
            excel_file %>%
            select(-all_of(columns_with_percents)) %>%
            pivot_longer(
              cols = all_of(columns_with_counts),
              names_to = "Marker_Full", values_to = "Count"
            ) %>%
            mutate(
              CellType = tab_names[tab],
              SourceFile = filepath#,
            ) %>%
            select(all_of(c(setdiff(names(.), "Count"), "Count"))) # *
          # * convoluted, but makes percent last column for tidiness
        }
      ) %>%
      bind_rows()
    
    if (!is.na(metadata_tab_num)) {
      pid <- read_excel(filepath, sheet = metadata_tab_num)
      return(list(Percents = percents, Counts = counts, Metadata = pid))
    } else {
      return(list(Percents = percents, Counts = counts))
    }
  }


# extract_markers_from_multiple_xlsx()
# ------------------------------------------------------------------------------
# if data is across multiple .xlsx sheets, automatically join
# input: multiple file paths as a string vector
# output: list as above
# ------------------------------------------------------------------------------
extract_markers_from_multiple_xlsx <-
  function(multiple_filepaths, metadata_tab_num = NULL) {
    if (is.null(metadata_tab_num)) {
      output_premerged <-
        map(multiple_filepaths,
          extract_flowcy_from_single_xlsx,
          metadata_tab_num = metadata_tab_num
        )
    } else {
      output_premerged <-
        tibble(
          filepath = multiple_filepaths,
          metadata_tab_num = metadata_tab_num
        ) %>%
        pmap(extract_flowcy_from_single_xlsx)
    }

    list(
      Percents = output_premerged %>% map_dfr(~ .$Percents),
      Counts = output_premerged %>% map_dfr(~ .$Counts),
      Metadata = output_premerged %>% map_dfr(~ .$Metadata)
    )
  }




# functions for Table 1 / plots ================================================
# (look excessive b/c copied from existing personal utility package 'choomisc')

# continuous variables
summaryMedianRange <-
  function(x, digits = 1, na.rm = F,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = T,
           symbol = NULL) {
    if (!(typeof(x) %in% c("double", "integer"))) {
      warning("input vector not numeric")
      x <- as.numeric(x)
    }


    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100

    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }

    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }

    output <-
      c(
        min(x, na.rm = na.rm),
        median(x, na.rm = na.rm),
        max(x, na.rm = na.rm)
      )

    output <- formatC(output, format = "f", digits = digits)

    output <- ifelse(is.null(symbol),
      paste0(output[2], " (", output[1], ", ", output[3], ")"),
      paste0(output[2], symbol, " (", output[1], symbol, ", ", output[3], symbol, ")")
    )

    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }

    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }

    return(output)
  }

# categorical variables
summaryCountPercent <-
  function(x, values, count_NAs = F,
           digits = 1,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = T) {
    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100

    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }

    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }

    n <- sum(x %in% values, na.rm = T)
    tot <- ifelse(count_NAs, length(x), length(x) - n_NA)

    output <- paste0(
      n, " (", formatC(n / tot * 100, format = "f", digits = digits),
      "%)"
    ) 

    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }

    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }

    return(output)
  }

# por plotting IQR
fxnMakeQuantiles <-
  function(x) {
    tmp <- quantile(x, c(0.25, 0.5, 0.75))
    tmp <- as.data.frame(t(tmp))
    names(tmp) <- c("ymin", "y", "ymax")
    tmp
  }



