#' Identify significant covariates with FWER, FDR, or no multiple testing
#' adjustment in the survival-data setting.
#'
#' Calculate p-values using pooled or per variable approach and identify
#' significant variables using FWER, FDR or no multiple testing adjustment of
#' p-values for a given alpha level of significance. This is the survival
#' equivalent of `add_test_results()`.
#'
#' @param vimpermsim List, an output of the `vim_perm_sim_survival()` function
#'   that stores `niters` variable importance values for both the original and
#'   row-wise permuted predictors.
#' @param alpha Numeric, the significance level, must be between 0 and 1.
#' @param init_num_vars Numeric, the number of covariates originally included in
#'   the data. Required to correctly apply the Benjamini-Hochberg (FDR) and Holm
#'   (FWER) p-value correction
#' @param to_show Character, one of `"FWER"`, `"FDR"` or `"unadjusted"`.
#'  * `"FWER"` (the default) - the output includes unadjusted,
#'   Benjamini-Hochberg (FDR) and Holm (FWER) adjusted p-values together with
#'   the decision whether the variable is significant or not (1 - significant, 0
#'   means not significant) according to the chosen criterium.
#'  * `"FDR"` - the output includes both unadjusted and FDR adjusted p-values
#'   along with the decision.
#'  * `"unadjusted:` - the output contains only raw, unadjusted p-values together
#'   with the decision.
#' @return A list with at least three elements. Its exact length depends on the
#' number of event levels in the input data (e.g., if status has levels
#' censoring, event, competing event, the list will have 4 elements: one data
#' frame per non-censoring event plus two lists).
#'  * Per-event data frames: For each non-censoring event in the input, the output
#'    contains a data frame named after that event. Example: if events are coded
#'    as 0 = censoring, 1 - event, and 2 - competing risk, the object includes two
#'    data frames: `event.1` and `event.2`. Each per-event data frame reports the
#'    variable and shadow importances produced by `vim_perm_sim_survival()`.
#'  * `test_res_pooled` - a list with one data frame per event, containing
#'    p-values computed using the pooled approach.
#'  * `test_res_per_variable` - a list with one data frame per event, containing
#'    p-values computed using the per-variable approach.
#'  All data frames within test_res_pooled and test_res_per_variable also include
#'  decision columns that flag variable importance based on the displayed
#'  p-values. Which decisions are shown (e.g., FWER, FDR, or unadjusted) is
#'  controlled by the to_show parameter.
#' @noRd
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust median ecdf sd
#' @examples
#' # Standard survival data: Veterans' Administration Lung Cancer study data
#' data(veteran, package = "randomForestSRC")
#' # Create vimpermsim object first
#' # When working with real data, increase num.trees value or leave default
#' # Here this parameter is set to a small value in order to reduce the run time
#' veteran_vps <- vim_perm_sim_survival(data = veteran,
#'  time_column = "time",
#'  status_column = "status",
#'  niters = 30,
#'  num.trees = 10)
#'
#' init_num_vars <- ncol(veteran) - 2
#'
#' # Display decisions based on all available p-values (FWER, FDR, unadjusted)
#' veteran_add_fwer <- add_test_results_survival(
#'   vimpermsim = veteran_vps, alpha = 0.05,
#'   init_num_vars = init_num_vars
#' )
#'
#' # Display decisions based on FDR adjusted and unadjusted p-values,
#' # expected warning
#' veteran_add_fdr <- suppressWarnings(add_test_results(
#'   vimpermsim = veteran_vps, alpha = 0.05,
#'   init_num_vars = init_num_vars, to_show = "FDR"
#' ))
#'
#' # Display decisions based on unadjusted p-values (Type1_confirmed column),
#' # expected warning
#' veteran_add_unadj <- suppressWarnings(add_test_results(
#'   vimpermsim = veteran_vps, alpha = 0.05,
#'   init_num_vars = init_num_vars, to_show = "unadjusted"
#' ))
add_test_results_survival <- function(vimpermsim,
                                      alpha = 0.05,
                                      init_num_vars,
                                      to_show = c("FWER", "FDR", "unadjusted")) {
  # Ensure to_show is one of the allowed values
  to_show <- match.arg(to_show)

  # Ensure alpha is numeric
  if (is.numeric(alpha) == FALSE) {
    stop("`alpha` must be numeric.")
  } else if (length(alpha) > 1 || alpha > 1 || alpha < 0) {
    stop("`alpha` must be a single numeric between 0 and 1.")
  }

  # Helper function
  # compare medians to pooled quantiles
  quants_pooled_fun <- function(varnames_originals, shadows, medians) {
    # compute standard deviations to make shadows of variables with different cardinalities comparable
    sd_of_shadows <- sapply(shadows, sd)
    mean_of_shadows <- sapply(shadows, mean)

    # divide medians and shadows by sd of shadow
    medians_divided_sd <- (medians - mean_of_shadows) / sd_of_shadows
    shadows_centered <- sweep(shadows, 2, mean_of_shadows, FUN = "-")
    shadows_divided_sd <- sweep(shadows_centered, 2, sd_of_shadows, FUN = "/")

    # ecdf function
    ec_pool <- ecdf(c(unlist(shadows_divided_sd), Inf))

    # return for all varnames remaining
    return(lapply(varnames_originals, function(x) {
      ec_pool(medians_divided_sd[[x]])
    }))
  }

  # Get the names of competing risks
  event_names <- names(vimpermsim)

  # Iterate over all competing risks
  for (event in event_names) {
    # Extract data for a selected competing risk
    vim_simulated <- vimpermsim[[event]]

    # split data
    originals <- vim_simulated %>%
      select(-ends_with("_permuted"))
    shadows <- vim_simulated %>%
      select(ends_with("_permuted"))

    # calculate medians of original variables
    medians <- originals %>% summarise(across(everything(), median))
    varnames_originals <- names(medians)

    # per variable
    # compare medians to quantiles of respective shadow variable
    quants_per_variable <- lapply(varnames_originals, function(x) {
      ec <- ecdf(c(shadows[[paste0(x, "_permuted")]], Inf))
      ec(medians[[x]])
    })

    # p-value correction - per variable method
    quants_per_variable_df <- data.frame(
      "varname" = varnames_originals,
      "quantile_per_variable" = unlist(quants_per_variable)
    ) %>%
      arrange(desc(.data[["quantile_per_variable"]])) %>%
      mutate(
        p_unadj = 1 - .data[["quantile_per_variable"]],
        p_adj_FDR = p.adjust(1 - .data[["quantile_per_variable"]], "BH", n = init_num_vars),
        p_adj_FWER = p.adjust(1 - .data[["quantile_per_variable"]], "holm", n = init_num_vars)
      ) %>%
      mutate(
        Type1_confirmed = ifelse(.data[["p_unadj"]] <= alpha, 1, 0),
        FDR_confirmed = ifelse(.data[["p_adj_FDR"]] <= alpha, 1, 0),
        FWER_confirmed = ifelse(.data[["p_adj_FWER"]] <= alpha, 1, 0)
      )

    # pooled
    quants_pooled <- quants_pooled_fun(varnames_originals, shadows, medians)

    # p-value correction - pooled method
    quants_pooled_df <- data.frame(
      "varname" = varnames_originals,
      "quantile_pooled" = unlist(quants_pooled)
    ) %>%
      arrange(desc(.data[["quantile_pooled"]])) %>%
      mutate(
        p_unadj = 1 - .data[["quantile_pooled"]],
        p_adj_FDR = stats::p.adjust(1 - .data[["quantile_pooled"]], "BH", n = init_num_vars),
        p_adj_FWER = stats::p.adjust(1 - .data[["quantile_pooled"]], "holm", n = init_num_vars)
      ) %>%
      mutate(
        Type1_confirmed = ifelse(.data[["p_unadj"]] <= alpha, 1, 0),
        FDR_confirmed = ifelse(.data[["p_adj_FDR"]] <= alpha, 1, 0),
        FWER_confirmed = ifelse(.data[["p_adj_FWER"]] <= alpha, 1, 0)
      )

    all_cols_names <- colnames(quants_pooled_df)

    if (to_show == "unadjusted") {
      to_select <- grep("varname|unadj|Type1", all_cols_names, value = TRUE)

      quants_pooled_df <- quants_pooled_df %>%
        select(all_of(c(to_select, "quantile_pooled")))

      quants_per_variable_df <- quants_per_variable_df %>%
        select(all_of(c(to_select, "quantile_per_variable")))

      warning("By setting the `to_show` parameter to `unadjusted` you will not be able to use the `plot_vimps()` function.")
    } else if (to_show == "FDR") {
      to_select <- grep("varname|unadj|Type1|FDR", all_cols_names, value = TRUE)

      quants_pooled_df <- quants_pooled_df %>%
        select(all_of(c(to_select, "quantile_pooled")))

      quants_per_variable_df <- quants_per_variable_df %>%
        select(all_of(c(to_select, "quantile_per_variable")))

      warning("By setting the `to_show` parameter to `FDR` you will not be able to use the `plot_vimps()` function.")
    }

    vimpermsim$test_res_pooled[[event]] <- quants_pooled_df
    vimpermsim$test_res_per_variable[[event]] <- quants_per_variable_df
  }

  return(vimpermsim)
}
