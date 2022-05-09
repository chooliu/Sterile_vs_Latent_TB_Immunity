# ==============================================================================
# 06_stimulated_comparisons.R
# same as 05, but on stim % positivity values
# ==============================================================================




# regression wrap functions ====================================================

# beta regression wrappers that return model objects
# similar motivation as in script 04, but now
# need +/- baseline adjustment functionality (Mtb-memory versus Mtb-stimulation)
# "nbl" = nobaseline

run_betareg_stim <-
  function(input_data) {
    betareg(Outcome ~ Group + Age + Percent_NoStim + Run,
            data = input_data, link = "logit")
  }

run_betareg_stim_norun <-
  function(input_data) {
    betareg(Outcome ~ Group + Age + Percent_NoStim,
            data = input_data, link = "logit")
  }

run_betareg_stim_nbl <-
  function(input_data) {
    betareg(Outcome ~ Group + Age + Run,
            data = input_data, link = "logit")
  }

run_betareg_stim_norun_nbl <-
  function(input_data) {
    betareg(Outcome ~ Group + Age,
            data = input_data, link = "logit")
  }

# return residuals
get_residuals_stim <-
  function(input_data) {
    betareg(Outcome ~ Age + Percent_NoStim + Run,
            data = input_data, link = "logit") %>%
      residuals(type = "deviance") %>%
      set_names(input_data$PID)
  }

get_residuals_stim_norun <-
  function(input_data) {
    betareg(Outcome ~ Age + Percent_NoStim,
            data = input_data, link = "logit") %>%
      residuals(type = "deviance") %>%
      set_names(input_data$PID)
  }

get_residuals_stim_nbl <-
  function(input_data) {
    betareg(Outcome ~ Age + Run,
            data = input_data, link = "logit") %>%
      residuals(type = "deviance") %>%
      set_names(input_data$PID)
  }

get_residuals_stim_norun_nbl <-
  function(input_data) {
    betareg(Outcome ~ Age,
            data = input_data, link = "logit") %>%
      residuals(type = "deviance") %>%
      set_names(input_data$PID)
  }




# master modeling function for stim ============================================
# similar arguments as in script 05 for unstim
# note new "correct_for_baseline" variable (i.e., 'subtract' unstim correction)

run_stim_analysis <-
  function(groups_to_include, group_model_coef,
           numerator, denominator,
           correct_for_run, invert_OR_dir, correct_for_baseline,
           inputdeltadat) {
    
    # filter stimulated data
    input_stim <-
      inputdeltadat %>%
      filter(Group %in% groups_to_include)
    
    # extract sample size for offsetting 0 and 100% values
    # (per betareg package documentation)
    sample_size <-
      input_stim %>%
      .$PID %>%
      unique %>%
      length
    
    input_stim %<>%
      mutate(Outcome = Percent_CMStim / 100,
             Outcome = (Outcome * (sample_size - 1) + 0.5) / sample_size) %>%
      group_by(Panel, CellType, Marker, FreqOf) %>%
      bind_cols(Group_Index = group_indices(.))
      
    # data driven filtering:
    # exclude markers with stim - unstim = |delta| <0.1% in >= 33% of subjects
    markersummary_stim <-
      input_stim %>%
      group_by(Panel, CellType, Marker) %>%
      summarize(percent_difference_pointonepercent = # [*]
                  sum(abs(Difference) < 0.1) / length(Percent_CMStim),
                percent_difference_onepercent = # 
                  sum(abs(Difference) < 1) / length(Percent_CMStim),
                percent_stim_onepercent = # [*]
                  sum(abs(Difference) < 1) / length(Percent_CMStim),
                IQR_Diff = quantile(abs(Difference), 0.75) -
                  quantile(abs(Difference), 0.25),
                IQR_CMStim = quantile(abs(Percent_CMStim), 0.75) -
                  quantile(abs(Percent_CMStim), 0.25),
                Group_Index = first(Group_Index),
                Exclude_Marker = any(Exclude_Marker),
                Exclude_CellType = any(Exclude_CellType)) %>%
      mutate(Include_DatDriven =
               percent_difference_pointonepercent < 0.33)

    valid_group_indices <-
      markersummary_stim %>%
      filter(Include_DatDriven & !Exclude_CellType & !Exclude_Marker) %>%
      .$Group_Index

    # apply filter to input data    
    input_stim <-
      input_stim %>%
      filter(Group_Index %in% valid_group_indices) %>%
      mutate(OutcomeForPCA = log(Outcome / (1 - Outcome))) %>%
      mutate(Group = fct_drop(Group))
    
    prediction_set_df <-
      prediction_set_df %>% filter(Group %in% groups_to_include)
    
    # no-run models (one is run in every case)
    if (correct_for_baseline) {
      stimulated_models_norun <-
        input_stim %>%
        group_split() %>%
        map(run_betareg_stim_norun)
      model_residuals <-
        input_stim %>%
        group_split() %>%
        map(get_residuals_stim_norun)
    } else {
      stimulated_models_norun <-
        input_stim %>%
        group_split() %>%
        map(run_betareg_stim_norun_nbl)
      model_residuals <-
        input_stim %>%
        group_split() %>%
        map(get_residuals_stim_norun_nbl)
    }
      
    # run-adjusted models (not run for BCG comparisons)
    # if run, overwrite other residual set
    if (correct_for_run) {
      if (correct_for_baseline) {
      stimulated_models <-
          input_stim %>%
          group_split() %>%
          map(run_betareg_stim)
      model_residuals <-
        input_stim %>%
        group_split() %>%
        map(get_residuals_stim)
      } else {
        stimulated_models <-
          input_stim %>%
          group_split() %>%
          map(run_betareg_stim_nbl)
        model_residuals <-
          input_stim %>%
          group_split() %>%
          map(get_residuals_stim_nbl)
      }
    }
    
    # obj to compile final results output
    stimulated_models_summary <-
      group_keys(input_stim)
    
    if (correct_for_run) {
      core_results_set <- stimulated_models
    } else {
      core_results_set <- stimulated_models_norun
    }
    
    # predicted values for hypothetical "average" sample
    # (more important in past iterations of data, where there were run effects)
    model_predictions <-
      core_results_set %>%
      map(function(x) {
        
        if (correct_for_run & correct_for_baseline) {
          unstim_vals_df <-
            x$model %>%
            group_by(Group, Run) %>%
            summarize(Percent_NoStim = mean(x$model$Percent_NoStim))
          prediction_set_df <-
            prediction_set_df %>%
            left_join(., unstim_vals_df)
        }
        
        if (!correct_for_run & correct_for_baseline) {
            unstim_vals_df <-
              x$model %>%
              group_by(Group) %>%
              summarize(Percent_NoStim = mean(x$model$Percent_NoStim))
            prediction_set_df <-
              prediction_set_df %>%
              left_join(., unstim_vals_df)
        }
        
        # prediction_set_df <-
        #   x$model
        
        # effect size & p-value for group
        predict_out <-
          tibble(Group = prediction_set_df$Group,
                 Predicted = predict(x, prediction_set_df)) %>%
          group_by(Group) %>%
          summarize(Predicted = mean(Predicted) * 100)
        
        predict_out$Predicted %>% set_names(predict_out$Group)
        
      }) %>%
      bind_rows() %>%
      set_names(., paste0("Predicted % ", names(.)))
    
    stimulated_models_summary$P_group <-
      core_results_set %>%
      map_dbl( function(x) {
        x %>% summary %>% .$coefficients %>% .$mean %>% .[group_model_coef, 4]
      })
    
    stimulated_models_summary$FDR_group <-
      p.adjust(stimulated_models_summary$P_group, "fdr")
    
    stimulated_models_summary$beta_group <-
      core_results_set %>%
      map_dbl( function(x) {
        x %>% summary %>% .$coefficients %>% .$mean %>% .[group_model_coef, 1]
      })
    
    stimulated_models_summary$OR_group <-
      exp(stimulated_models_summary$beta_group)
    
    stimulated_models_summary$OR_CI_low <- 
      core_results_set %>%
      map_dbl( function(x) {
        confint(x) %>%
          .[2, "2.5 %"] %>%
          exp })
    stimulated_models_summary$OR_CI_high <- 
      core_results_set %>%
      map_dbl( function(x) {
        confint(x) %>%
          .[2, "97.5 %"] %>%
          exp })
    
    if (invert_OR_dir) {
      stimulated_models_summary <-
        stimulated_models_summary %>%
        mutate_at(.vars = c("OR_group", "OR_CI_low", "OR_CI_high"),
                  .funs = function(x) { 1/x } )
    }
    
    
    stimulated_models_summary$P_age <-
      core_results_set %>%
      map_dbl( function(x) {
        x %>% summary %>% .$coefficients %>% .$mean %>%
          .["Age", 4]
      })
    
    if (correct_for_run) {
      stimulated_models_summary$P_run <-
        seq(1:nrow(stimulated_models_summary)) %>%
        map_dbl( function(i) {
          waldtest(stimulated_models[[i]], stimulated_models_norun[[i]]) %>%
            .$`Pr(>Chisq)` %>% .[2]
        } )
    } else {
      stimulated_models_summary$P_run <- NA
    }
    
    # final results summary
    order_of_pvals <- order(stimulated_models_summary$P_group)
    
    resultsout_pt1 <-
      stimulated_models_summary %>%
      rowwise() %>%
      transmute(Panel, CellType, Marker,
                `Signif (FDR < 0.05)` = if_else(FDR_group < 0.05, "*", "")) %>%
      bind_cols(model_predictions)
    
    resultsout_pt2 <-
      stimulated_models_summary %>%
      rowwise() %>%
      transmute(`Raw Stim % in Numerator Group` =
                  calc_raw_dat_percent(Panel, CellType, Marker, numerator,
                                       "Stim", inputdeltadat),
                `Raw Stim % in Denominator Group` =
                  calc_raw_dat_percent(Panel, CellType, Marker, denominator,
                                       "Stim", inputdeltadat),
                `Raw Delta % in Numerator Group` =
                  calc_raw_dat_percent(Panel, CellType, Marker, numerator,
                                       "Delta", inputdeltadat),
                `Raw Delta % in Denominator Group` =
                  calc_raw_dat_percent(Panel, CellType, Marker, denominator,
                                       "Delta", inputdeltadat),
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
      unite(col = "OR (95% CI)", `OR (2.5%)`, `OR (97.5%)`,
            sep = " - ", remove = T) %>%
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
         Residuals = model_residuals %>% bind_cols(input_stim %>% group_keys, .),
         Predicted = model_predictions %>% bind_cols(input_stim %>% group_keys, .)) %>%
      map(~ .x[order_of_pvals, ])
  }




# run baseline-adjusted models, functionality ==================================
# "Mtb-memory" results in manuscript

resultstable_stim_lneg_lpos <-
  run_stim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                    correct_for_run = F, invert_OR_dir = F, correct_for_baseline = T, 
                    deltas_final %>% filter(FreqOf == "Parent"))
resultstable_stim_lneg_bcg <-
   run_stim_analysis(c("LTBI-", "BCG"), "GroupLTBI-", "BCG", "LTBI-",
                     correct_for_run = F, invert_OR_dir = T, correct_for_baseline = T, 
                     deltas_final %>% filter(FreqOf == "Parent"))
resultstable_stim_lpos_bcg <-
  run_stim_analysis(c("LTBI+", "BCG"), "GroupLTBI+", "BCG", "LTBI+",
                    correct_for_run = F, invert_OR_dir = T, correct_for_baseline = T, 
                    deltas_final %>% filter(FreqOf == "Parent"))




# run no baseline-adjusted models ==============================================
# (was formerly run as a sensitivity analysis,
# now excluded from final paper because creates confusion)

resultstable_stim_lneg_lpos_nbl <-
  run_stim_analysis(c("LTBI-", "LTBI+"), "GroupLTBI+", "LTBI+", "LTBI-",
                    correct_for_run = F, invert_OR_dir = F, correct_for_baseline = F,
                    deltas_final %>% filter(FreqOf == "Parent"))
resultstable_stim_lneg_bcg_nbl <-
  run_stim_analysis(c("LTBI-", "BCG"), "GroupLTBI-", "BCG", "LTBI-",
                    correct_for_run = F, invert_OR_dir = T, correct_for_baseline = F,
                    deltas_final %>% filter(FreqOf == "Parent"))
resultstable_stim_lpos_bcg_nbl <-
  run_stim_analysis(c("LTBI+", "BCG"), "GroupLTBI+", "BCG", "LTBI+",
                    correct_for_run = F, invert_OR_dir = T, correct_for_baseline = F,
                    deltas_final %>% filter(FreqOf == "Parent"))




# summarize # tests, significant results =======================================

mget(ls() %>% .[grepl("resultstable_stim_", .)]) %>%
  map( ~ .x$Results %>% nrow() ) %>% unlist %>% sum

mget(ls() %>% .[grepl("resultstable_stim_", .)]) %>%
  map( ~ .x$Results$`Signif (FDR < 0.05)` %>% table() )



