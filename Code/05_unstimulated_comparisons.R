# ==============================================================================
# 05_unstimulated_comparisons.R
# use beta regression to compare unstimulated % positivity values between groups
# ==============================================================================




# regression wrap functions ====================================================

# beta regression wrappers that return model objects
# note that Group has to be first on RHS after ~ in model
# (some subsequent fxns are hard-coded to extract second coefficient)
# cauchy link function used since interest in markers close to 0% or 100%

run_betareg_unstim <-
  function(input_data) {
    print(unique(input_data$CellType))
    print(unique(input_data$Marker))
    betareg(Outcome ~ Group + Age + Run, data = input_data,
            control = betareg.control(maxit  = 10000),
            link = "logit")
  }

run_betareg_unstim_norun <-
  function(input_data) {
    print(unique(input_data$CellType))
    print(unique(input_data$Marker))
    betareg(Outcome ~ Group + Age, data = input_data,
            link = "logit")
  }

# return residuals
get_residuals_unstim <-
  function(input_data) {
    print(unique(input_data$CellType))
    print(unique(input_data$Marker))
    betareg(Outcome ~ Age + Run, data = input_data,
            link = "logit") %>%
      residuals(type = "deviance") %>%
      set_names(input_data$PID)
  }

get_residuals_unstim_norun <-
  function(input_data) {
    print(unique(input_data$CellType))
    print(unique(input_data$Marker))
    betareg(Outcome ~ Age, data = input_data,
            link = "logit") %>%
      residuals(type = "deviance") %>%
      set_names(input_data$PID)
  }

# return predicted: data for
# approximate group-effect assuming median age = 30
prediction_set_df <-
  tribble(~Group, ~Age, ~Run,
          "LTBI-", 7, "Run1",
          "LTBI+", 7, "Run1",
          "LTBI-", 7, "Run2",
          "LTBI+", 7, "Run2",
          "BCG", 7, NA_character_) %>%
  mutate(Group = as.factor(Group))

# for results reporting, report raw dat % values
calc_raw_dat_percent <-
  function(p_in, ct_in, m_in, group, stim, inputdeltadat) {
    
    inputdeltadat %>%
      filter(Panel == p_in & CellType == ct_in &
               Marker == m_in & Group == group) %>%
      summarize(NoStim = mean(Percent_NoStim),
                Stim = mean(Percent_CMStim),
                Delta = mean(Percent_CMStim - Percent_NoStim)) %>%
      .[stim] %>%
      unlist
  }




# master modeling function for unstim ==========================================
# - groups_to_include  (pairwise between-grp comparisons)
# - group_model_coef (e.g., "GroupLTBI+") [required bc helps with error checking]
# - numerator / denominator (for pretty outputs, clarify effect direction)
# - correct_for_run (add indicator variable for run)
# - invert_OR_dir (after analysis was performed, decided to make ref group BCG)

run_unstim_analysis <-
  function(groups_to_include, group_model_coef,
           numerator, denominator, correct_for_run,
           invert_OR_dir,
           inputdeltadat) {

    # filters unstimulated data
    input_unstim <-
      inputdeltadat %>%
      filter(Group %in% groups_to_include)

    # extract sample size for offsetting 0 and 100% values
    # (per betareg package documentation)
    sample_size <-
      input_unstim %>%
      .$PID %>%
      unique %>%
      length
    
    input_unstim %<>%
      mutate(Outcome = Percent_NoStim / 100,
             Outcome = (Outcome * (sample_size - 1) + 0.5) / sample_size) %>%
      group_by(Panel, CellType, Marker, FreqOf) %>%
      bind_cols(Group_Index = group_indices(.))
    
    # data driven filtering: exclude markers with <0% in >= 33% of subjects
    markersummary_unstim <-
      input_unstim %>%
      group_by(Panel, CellType, Marker) %>%
      summarize(percent_unstim_zeroes = # [*]
                  sum(Percent_NoStim == 0) / length(Percent_NoStim),
                Group_Index = first(Group_Index),
                Exclude_Marker = any(Exclude_Marker),
                Exclude_CellType = any(Exclude_CellType)) %>%
      mutate(Include_DatDriven =
               percent_unstim_zeroes < 0.33)
    
    valid_group_indices <-
      markersummary_unstim %>%
         filter(Include_DatDriven & !Exclude_CellType & !Exclude_Marker) %>%
      .$Group_Index
    
    # apply filter to input data
    input_unstim <-
      input_unstim %>%
      filter(Group_Index %in% valid_group_indices) %>%
      mutate(Group = fct_drop(Group))
    
    prediction_set_df <-
      prediction_set_df %>% filter(Group %in% groups_to_include) %>%
      filter(!duplicated(Group))
    
    # run models (beta regression)
    unstimulated_models_norun <-
      input_unstim %>%
      group_split() %>%
      map(run_betareg_unstim_norun)
    
    model_residuals <-
      input_unstim %>%
      group_split() %>%
      map(get_residuals_unstim_norun)
    
    if (correct_for_run) {
      unstimulated_models <-
        input_unstim %>%
        group_split() %>%
        map(run_betareg_unstim)
      model_residuals <-
        input_unstim %>%
        group_split() %>%
        map(get_residuals_unstim)
    }

    # obj to compile final results output
    unstimulated_models_summary <-
      group_keys(input_unstim)
    
    if (correct_for_run) {
      core_results_set <- unstimulated_models
    } else {
      core_results_set <- unstimulated_models_norun
    }
    
    # predicted values for hypothetical "average" sample
    # (more important in past iterations of data, where there were run effects)
    model_predictions <-
      core_results_set %>%
      map(function(x) {
        predict_out <-
          tibble(Group = prediction_set_df$Group,
               Predicted =
                 predict(x,
                         prediction_set_df)) %>%
          group_by(Group) %>%
          summarize(Predicted = mean(Predicted)*100)
          predict_out$Predicted %>% set_names(predict_out$Group)
      }) %>%
      bind_rows %>%
      set_names(., paste0("Predicted % ", names(.)))
    
    # effect size & p-value for group
    unstimulated_models_summary$P_group <-
      core_results_set %>%
      map_dbl( function(x) {
        x %>% summary %>% .$coefficients %>% .$mean %>% .[group_model_coef, 4]
      })
    
    unstimulated_models_summary$FDR_group <-
      p.adjust(unstimulated_models_summary$P_group, "fdr")
    
    unstimulated_models_summary$beta_group <-
      core_results_set %>%
      map_dbl( function(x) {
        x %>% summary %>% .$coefficients %>% .$mean %>% .[group_model_coef, 1]
      })
    
    unstimulated_models_summary$OR_group <-
      exp(unstimulated_models_summary$beta_group)
    unstimulated_models_summary$OR_CI_low <- 
      core_results_set %>%
      map_dbl( function(x) {
        confint(x) %>%
          .[2, "2.5 %"] %>%
          exp })
    unstimulated_models_summary$OR_CI_high <- 
      core_results_set %>%
      map_dbl( function(x) {
      confint(x) %>%
      .[2, "97.5 %"] %>%
      exp })
    
    if (invert_OR_dir) {
      unstimulated_models_summary <-
        unstimulated_models_summary %>%
        mutate_at(.vars = c("OR_group", "OR_CI_low", "OR_CI_high"),
                  .funs = function(x) { 1/x } )
    }
    
    # age and run effects
    unstimulated_models_summary$P_age <-
      core_results_set %>%
      map_dbl( function(x) {
        x %>% summary %>% .$coefficients %>% .$mean %>%
          .["Age", 4]
      })
    
    if (correct_for_run) {
    unstimulated_models_summary$P_run <-
      seq(1:nrow(unstimulated_models_summary)) %>%
      map_dbl( function(i) {
        waldtest(unstimulated_models[[i]], unstimulated_models_norun[[i]]) %>%
          .$`Pr(>Chisq)` %>% .[2]
      } )
    } else {
      unstimulated_models_summary$P_run <- NA
    }
  
    # final results summary
    order_of_pvals <- order(unstimulated_models_summary$P_group)
    
    resultsout_pt1 <-
      unstimulated_models_summary %>%
      rowwise() %>%
      transmute(Panel, CellType, Marker,
                `Signif (FDR < 0.05)` = if_else(FDR_group < 0.05, "*", "")) %>%
      bind_cols(model_predictions)
    
    resultsout_pt2 <-
      unstimulated_models_summary %>%
      rowwise() %>%
      transmute(`Raw Unstim % in Numerator Group` =
                  calc_raw_dat_percent(Panel, CellType, Marker, numerator, "NoStim", inputdeltadat),
                `Raw Unstim % in Denominator Group` =
                  calc_raw_dat_percent(Panel, CellType, Marker, denominator, "NoStim", inputdeltadat),
                `Odds Ratio (Odds Numerator / Odds Denominator)` = OR_group,
                `OR (2.5%)` = OR_CI_low,
                `OR (97.5%)` = OR_CI_high,
                `P (Numerator vs Denominator)` = P_group,
                `FDR (Numerator vs Denominator)` = FDR_group,
                `P (Age)` = P_age,
                `P (Run)` = P_run) %>%
      mutate_if(is.numeric, format, format = "fg", digits = 2, nsmall = 2) %>%
      mutate(`OR (2.5%)` = paste0("(", `OR (2.5%)`),
             `OR (97.5%)` = paste0(`OR (97.5%)`, ")")) %>%
      unite(col = "OR (95% CI)", `OR (2.5%)`, `OR (97.5%)`, sep = " - ", remove = T) %>%
      unite(col = "Odds Ratio (Odds Numerator / Odds Denominator)",
            `Odds Ratio (Odds Numerator / Odds Denominator)`, `OR (95% CI)`,
            sep = "\n", remove = T)
    
    resultsout <-
      bind_cols(resultsout_pt1, resultsout_pt2) %>%
      mutate_if(is.numeric, format, format = "fg", digits = 2, nsmall = 2)
    
    if (invert_OR_dir) {
      names(resultsout) %<>%
        gsub("Numerator", denominator, .) %>%
        gsub("Denominator", numerator, .)
    } else {
      names(resultsout) %<>%
        gsub("Numerator", numerator, .) %>%
        gsub("Denominator", denominator, .)
    }
    
    # export model residuals for visualization
    model_residuals <-
      model_residuals %>%
      do.call(bind_rows, .)
    
    list(Results = resultsout,
         Residuals = model_residuals %>% bind_cols(input_unstim %>% group_keys, .),
         Predicted = model_predictions %>% bind_cols(input_unstim %>% group_keys, .)) %>%
      map(~ .x[order_of_pvals, ])
  }




# run models for fxn % positivity ==============================================

resultstable_unstim_lneg_lpos <-
  run_unstim_analysis(c("LTBI+", "LTBI-"), "GroupLTBI+", "LTBI+", "LTBI-",
                      correct_for_run = F, invert_OR_dir = F,
                      inputdeltadat = deltas_final %>% filter(FreqOf == "Parent"))
resultstable_unstim_lneg_bcg <-
  run_unstim_analysis(c("LTBI-", "BCG"), "GroupLTBI-", "BCG", "LTBI-",
                      correct_for_run = F, invert_OR_dir = T,
                      inputdeltadat = deltas_final %>% filter(FreqOf == "Parent"))
resultstable_unstim_lpos_bcg <-
  run_unstim_analysis(c("LTBI+", "BCG"), "GroupLTBI+", "BCG", "LTBI+",
                      correct_for_run = F, invert_OR_dir = T,
                      inputdeltadat = deltas_final %>% filter(FreqOf == "Parent"))




# run models for population %s =================================================

resultstable_unstim_lneg_lpos_parent <-
  run_unstim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                      correct_for_run = F, invert_OR_dir = F,
                      deltas_final %>% filter(FreqOf != "Parent"))
resultstable_unstim_lneg_bcg_parent <-
  run_unstim_analysis(c("LTBI-", "BCG"), "GroupLTBI-", "BCG", "LTBI-",
                      correct_for_run = F, invert_OR_dir = T,
                      deltas_final %>% filter(FreqOf != "Parent"))
resultstable_unstim_lpos_bcg_parent <-
  run_unstim_analysis(c("LTBI+", "BCG"), "GroupLTBI+", "BCG", "LTBI+",
                      correct_for_run = F, invert_OR_dir = T,
                      deltas_final %>% filter(FreqOf != "Parent"))




# summarize # tests, significant results =======================================

mget(ls() %>% .[grepl("resultstable_unstim_", .)]) %>%
  map( ~ .x$Results %>% nrow() ) %>% unlist %>% sum

mget(ls() %>% .[grepl("resultstable_unstim_", .)]) %>%
  map( ~ .x$Results$`Signif (FDR < 0.05)` %>% table() )

