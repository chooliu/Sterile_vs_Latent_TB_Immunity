# ==============================================================================
# 01_load_metadata_and_markers.R
# load and clean metadata + activation marker .xlsx files &
# also create various long -> wide versions
# ==============================================================================




# load activation markers ======================================================

# percent positivity in wide format from E.J. flowjo export,
# samples (rows) x fraction of parent or raw event counts (columns)
exceldat <-
  extract_markers_from_multiple_xlsx(
    c(
      "./DataRaw/20210528/20190716 APC TB Orig unmixed v20210519 for stats.xlsx",
      "./DataRaw/20210819/20190716 TB Tcell v20210407 for stats v4.xlsx"
    ),
    c(10, 12)
  )

dat_frequencies <-
  exceldat$Percents %>%
  transmute(
    Sample = toupper(`Sample:`),
    PID = if_else(is.na(PID),
      str_split(Sample, " ") %>% map_chr(~ .[2]), PID
    ),
    Panel = if_else(
      grepl("APC", SourceFile), "APC", "TCell"
    ),
    CellType, Marker_Full, Percent
  ) %>%
  filter(!is.na(Percent) & !is.na(PID))

dat_counts <-
  exceldat$Counts %>%
  transmute(
    Sample = toupper(`Sample:`),
    PID = if_else(is.na(PID),
      str_split(Sample, " ") %>% map_chr(~ .[2]), PID
    ),
    Panel = if_else(
      grepl("APC", SourceFile), "APC", "TCell"
    ),
    CellType, Marker_Full, Count
  ) %>%
  filter(!is.na(Count) & !is.na(PID))




# format marker names ==========================================================
# goal is to go from long gating tree name (via FlowJo, long strings)
# shorten 'marker' name into a non-redundant clean celltype & marker
# e.g., "Singlets/Live/Size/Cleanup/dump-/DR+/CD14+/CD40" --> CD40s in CD14+

generate_shortname <-
  function(markerfull, add_end_depth = 0) {
    markerfull %>%
      str_split(., pattern = "/") %>%
      map_chr(.f = function(string) {
        stringleng <- length(string)
        paste0(string[seq(stringleng - add_end_depth, stringleng)],
          collapse = "/"
        )
      })
  }

# check if a given depth (# of "/"s from right) makes markers non-redundant
# if export_table = F, just checks yes/no redundancy,
# if export_table = T, prints putative celltype/marker combos
check_shortname_consistency <-
  function(freq_df, add_end_depth = 0, export_table = F) {
    outtable <-
      freq_df %>%
      mutate(
        Marker =
          str_split(Marker_Full, pattern = "/") %>%
            map_chr(.f = function(string) {
              stringleng <- length(string)
              paste0(string[seq(stringleng - add_end_depth, stringleng)],
                collapse = "/"
              )
            })
      ) %>%
      group_by(CellType)

    if (!export_table) {
      outtable <-
        outtable %>%
        summarize_at(
          .vars = c("Marker_Full", "Marker"),
          .funs = function(x) {
            length(unique(x))
          }
        ) %>%
        filter(Marker_Full != Marker | is.na(Marker))

      if (nrow(outtable) == 0) {
        print("non-redundant depth")
      } else {
        print("some redundant marker names")
        return(outtable)
      }
    } else {
      outtable %>%
        summarize_at(
          .vars = c("Marker"),
          .funs = list(
            `# Markers` = function(x) {
              length(unique(x))
            },
            `Unique Markers` = function(x) {
              unique(x) %>%
                gsub(" \\| Freq. of Parent| \\| Count| \\| Freq. of Cleanup| \\| Freq. of DR\\+", "", .) %>%
                paste0(collapse = ";   ")
            }
          )
        )
    }
  }


# the "freq of live (cleanup)" cells have some markers that need depth 2
dat_frequencies %>% check_shortname_consistency()
dat_frequencies %>% check_shortname_consistency(add_end_depth = 1)
dat_counts %>% check_shortname_consistency()

dat_frequencies <-
  dat_frequencies %>%
  mutate(Marker = if_else(CellType == "Parent Pop freq of DR" |
    CellType == "Cell Pop freq of live (cleanup)",
  generate_shortname(Marker_Full, 1),
  generate_shortname(Marker_Full, 0)
  )) %>%
  separate(col = Marker, into = c("Marker", "FreqOf"), sep = " \\| Freq. of ")

dat_counts <-
  dat_counts %>%
  mutate(Marker = generate_shortname(Marker_Full, 0) %>%
    gsub(" \\| Count", "", .))




# merge sample metadata ========================================================

# metadata of each sample compiled from 5 different spreadsheets by me;
# script generating the below .xlsx hidden to avoid showing participant metadata
metadata_consolated <-
  read_excel("./DataRaw/20210602/20210601_TB_compiled_metadat_CL EJ.xlsx") %>%
  dplyr::select(-Panel, -Sample, -Stim) %>%
  mutate(Group = as.factor(Group), Age = as.numeric(Age))


# check that marker IDs can be merged with metadata IDs:
# are all samples in flowho data in metadata & vice versa?
setdiff(dat_frequencies$PID, metadata_consolated$PID) %>% unique() # none
setdiff(metadata_consolated$PID, dat_frequencies$PID) %>% unique() # some BCGs





# finaliize data preprocessing =================================================
# (i) join marker & participant metadata (joined_data & frequencies_final)
# (ii) final manuscript names & subset examination
# (iii) make wide versions (deltas_final) for easier STIM - UNSTIM calcs

# (i) joined data in long format (22088 x 26) ----------------------------------
# * infer stim / unstim from sample name
joined_data <-
  dat_frequencies %>%
  left_join(., metadata_consolated, by = c("PID")) %>%
  mutate(Stim = if_else(grepl("CM", Sample), "CMStim", "NoStim")) # *


# (ii) final MS names ----------------------------------------------------------
# export spreadsheet for collaborator review:
# make sure greek symbols appear in manuscript-ready fmt (e.g., "IFN-gamma"),
# common language like "CD1c DC" --> cDC2,
# and column for a priori exclusion of subsets

# # export processed celltype/marker combos from line 43 for AW and EJ to review
# tibble(FlowJo =
#          c(dat_frequencies$CellType, dat_frequencies$Marker,
#            dat_counts$CellType, dat_counts$Marker) %>%
#          unique) %>%
#   mutate(Manuscript = "", Exclude = "") %>% # collaborators' MS name/exclusion
#   write_xlsx(path = "./Reports/20210922_TB_Manuscript_Names.xlsx",
#              format_headers = F)

# import filled out version of above containing final names
tidy_ms_names <-
  read_excel("DataRaw/20210819/20210922_TB_Manuscript_Names_CL_AW_EJ.xlsx") %>%
  transmute(FlowJo,
    Manuscript = if_else(is.na(Manuscript), FlowJo, Manuscript),
    Exclude, Notes = `...4`
  ) %>%
  # few manual changes: fix il1b specifically, exclude dnDCs
  mutate(
    Manuscript = if_else(Manuscript == "IL1b", "IL1β", Manuscript),
    Exclude = if_else(Manuscript == "cDCn", T, Exclude)
  )

# join with MS-ready names
frequencies_final <-
  joined_data %>%
  left_join(tidy_ms_names %>%
    rename(CellType_MS = Manuscript, Exclude_CellType = Exclude),
  by = c("CellType" = "FlowJo")
  ) %>%
  left_join(tidy_ms_names %>%
    rename(Marker_MS = Manuscript, Exclude_Marker = Exclude),
  by = c("Marker" = "FlowJo")
  )

# (ii) final MS names ----------------------------------------------------------

# # make celltypes & marker names publication-ready
# # if worked, the following should return nothing (i.e., no missing MS names)
frequencies_final %>%
  filter(is.na(Marker_MS) | is.na(CellType_MS)) %>%
  group_by(Marker, CellType) %>%
  group_keys()

frequencies_final <-
  frequencies_final %>%
  mutate(CellType = CellType_MS, Marker = Marker_MS) %>%
  filter(Group %in% c("BCG", "LTBI-", "LTBI+")) %>%
  mutate(Group = fct_drop(Group))

# "deltas_final" primary tibble used throughout; wide fmt
# note: in final manuscript LTBI- participants called "TB-resisters"
deltas_final <-
  frequencies_final %>%
  pivot_wider(
    id_cols = c(
      Panel, PID, Marker, CellType, FreqOf,
      Exclude_Marker, Exclude_CellType
    ),
    names_from = Stim,
    names_prefix = "Percent_",
    values_from = Percent
  ) %>%
  left_join(metadata_consolated, by = "PID") %>%
  filter(Run != "Run4") %>% # make sure run 4 is excluded; detected batch effects
  mutate(Difference = Percent_CMStim - Percent_NoStim) %>% # delta mainly for plotting
  filter(Group %in% c("BCG", "LTBI-", "LTBI+")) %>% # remove samples marked "exclude"
  mutate(
    Group = fct_drop(Group),
    GroupShort = fct_recode(Group, `L-` = "LTBI-", `L+` = "LTBI+", B = "BCG"),
    RunShort = gsub("Run", "", Run)
  ) # short versions for plotting




# export marker summary ========================================================

# (i) summarize raw event counts (though not analyzed; focus on frequencies)
# (ii) final list of celltypes/markers after renaming/preprocessing

parent_count_summary <-
  dat_counts %>%
  mutate(Stim = if_else(grepl("CM", Sample), "CMStim", "NoStim")) %>%
  group_by(Panel, CellType) %>%
  left_join(tidy_ms_names %>%
    rename(
      CellType_MS = Manuscript,
      Exclude_CellType = Exclude
    ),
  by = c("CellType" = "FlowJo")
  ) %>%
  mutate(CellType = CellType_MS) %>%
  left_join(metadata_consolated, by = "PID") %>%
  filter(PID %in% unique(deltas_final$PID)) %>%
  mutate(Group = fct_drop(Group)) %>%
  pivot_wider(
    id_cols = c(Panel, PID, CellType, Group),
    names_from = Stim,
    names_prefix = "Count_",
    values_from = Count
  ) %>%
  group_by(CellType, Group) %>%
  summarize(
    `Unstim%Less30` = sum(Count_NoStim < 30) / n(),
    `UnstimMin` = min(Count_NoStim),
    `UnstimQ25` = quantile(Count_NoStim, 0.25),
    `UnstimMedian` = quantile(Count_NoStim, 0.5),
    `UnstimQ75` = quantile(Count_NoStim, 0.75),
    `UnstimMax` = max(Count_NoStim),
    `Stim%Less30` = sum(Count_CMStim < 30) / n(),
    `StimMin` = min(Count_CMStim),
    `StimQ25` = quantile(Count_CMStim, 0.25),
    `StimMedian` = quantile(Count_CMStim, 0.5),
    `StimQ75` = quantile(Count_CMStim, 0.75),
    `StimMax` = max(Count_CMStim),
  ) %>%
  pivot_longer(cols = 3:ncol(.)) %>%
  pivot_wider(
    id_cols = c("CellType"),
    names_from = c("name", "Group"), names_sep = " "
  )

parsed_features_freq <-
  deltas_final %>%
  mutate(
    Marker = paste0(Marker, if_else(Exclude_Marker, "[EX]", "")),
    CellType = paste0(CellType, if_else(Exclude_CellType, "[EX]", ""))
  ) %>%
  group_by(CellType) %>%
  summarize(
    `# Markers` = unique(Marker) %>% .[!grepl("\\[EX", .)] %>% length(),
    `Unique Markers` = unique(Marker) %>% paste0(collapse = ";   ")
  ) %>%
  arrange(!grepl("APC|PBMC", CellType))



