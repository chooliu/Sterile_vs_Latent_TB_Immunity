# ==============================================================================
# 10_SPICE_fdr_correct.R
# FDR-adjustment for permutation p-values from SPICE
# ==============================================================================




# import SPICE-values ==========================================================
# i manually entered p-values from SPICE
# (p-vals sent to me in hard-coded figure form --> export p-values)

SPICE_pval_thresholds <-
  read_excel("./DataRaw/20220302/20220305_SPICE_transcription.xlsx") %>%
  apply(., 2,
        function(pvals) {
          
          # proper adjustment of <0.00001 p-value for number of permutations
          # for data entry entered <0.0001 values from SPICE as "0"
          # but change 0 --> 1/10001, which is from the # of permutations 10K + 1
          
          pvals <- pvals[!is.na(pvals)]
          pvals[pvals == 0] <- 1/10001
          
          # what's the largest p-value that corresponds to the fdr < 0.05 threshold?
          tibble(p = pvals,
                 fdr = p.adjust(pvals, method = "fdr")) %>%
            arrange(-fdr) %>%
            filter(fdr <= 0.05) %>%
            .[1, "p", drop = T]
  })

SPICE_pval_thresholds




# alternatively, global threshold ==============================================
# p <= 0.0247, across all the polyfunctionality tests

global_threshold <-
  read_excel("./DataRaw/20220302/20220305_SPICE_transcription.xlsx") %>%
  unlist %>%
  .[!is.na(.)] %>%
  sort

global_threshold <-
  max(global_threshold[p.adjust(global_threshold, "fdr") <= 0.05])




