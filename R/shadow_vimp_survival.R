#' Select influential covariates using random forest with multiple-testing
#' control in the survival-data setting.
#'
#' `shadow_vimp_survival()` performs variable selection and determines whether
#' each covariate is influential based on unadjusted, FDR-adjusted, and
#' FWER-adjusted p-values.This is the surviva equivalent of `shadow_vimp()`.
#'
#' @details
#' **Input data**:
#' The `shadow_vimp_survival()` function uses the survival random forest
#' implementation from \insertCite{rfpackage}{shadowVIMP}. Permissible input data
#' types are numeric (real-valued), integer, factor, or logical.
#'
#' **Pre-selection**:
#' The `shadow_vimp_survival()` function by default performs variable selection
#' in multiple steps. Initially, it prunes the set of predictors using a relaxed
#' (higher) alpha threshold in a pre-selection stage. Variables that pass this
#' stage then undergo a final evaluation using the target (lower) alpha
#' threshold and more iterations. This stepwise approach distinguishes
#' informative from uninformative covariates based on their VIMPs and enhances
#' computational efficiency. The user can also perform variable selection in a
#' single step, without a pre-selection phase.
#'
#' **Missing values**:
#' Missing values in both outcomes and predictors can be imputed using the
#' algorithm introduced in \insertCite{missing_val_alg}{shadowVIMP}. For details
#'  of the imputation algorithm, see the original paper or the survival-data
#'  vignette of this package.
#' @references
#' \insertAllCited{}
#' @param alphas Numeric vector, significance level values for each step of the
#'   procedure, default `c(0.3, 0.10, 0.05)`.
#' @param niters  Numeric vector, number of permutations to be performed in each
#'   step of the procedure, default `c(30, 120, 1500)`.
#' @param data Input data frame.
#' @param time_column Character, name of the column containing the event time.
#' @param status_column Character, name of the column containing the status
#'  indicator. The current implementation of the function supports the
#'  competing-risks setting.
#' @param to_show Character, one of `"FWER"`, `"FDR"` or `"unadjusted"`.
#'  * `"FWER"` (the default) - the output includes unadjusted,
#'   Benjamini-Hochberg (FDR) and Holm (FWER) adjusted p-values together with
#'   the decision whether the variable is significant or not (1 - significant, 0
#'   means not significant) according to the chosen criterium.
#'  * `"FDR"` - the output includes both unadjusted and FDR adjusted p-values
#'  along with the decision.
#'  * `"unadjusted:` - the output contains only raw, unadjusted p-values
#'   together with the decision.
#' @param num.trees Numeric, number of trees. Passed to [randomForestSRC::rfsrc()],
#'   default is `max(2 * (ncol(data) - 1), 10000)`.
#' @param importance Character, the type of variable importance to be calculated
#'   for each variable. Argument passed to [randomForestSRC::rfsrc()], default is
#'   `permutation`.
#' @param save_vimp_history Character, specifies which variable importance
#'   measures to save. Possible values are:
#'  * `"all"` (the default) - save variable importance measures from all steps
#'   of the procedure (both the pre-selection phase and the final selection
#'   step).
#'  * `"last"` - save only the variable importance measures from the final
#'   step.
#'  * `"none"` - do not save any variable importance measures.
#' @param method Character, one of `"pooled"` or `"per_variable"`.
#'  * `"pooled"` (the default) - the results of the final step of the procedure
#'  show the p-values obtained using the "pooled" approach and the corresponding
#'  decisions.
#'  * `"per_variable"` - the results of the final step of the procedure
#'  show the p-values obtained using the "per variable" approach and the
#'  corresponding decisions.
#' @param na.action Character, one of `"na.omit"` (default) or `"na.impute"`.
#' Action to take when the data contain `NA`. By default,
#' `shadow_vimp_survival()` excludes rows with missing values in any covariate.
#' When `na.action = "na.impute"`, missing values are imputed using the
#' imputation algorithm. See [randomForestSRC::rfsrc()] and the vignette on this
#' package’s survival-data extension for details.
#' @param ... Additional parameters passed to [randomForestSRC::rfsrc()].
#' @return Object of the class "shadow_vimp_survival" with the following entries:
#'  * `call` - the call formula used to generate the output.
#'  * `alpha` - numeric, significance level used in the algorithm.
#'  * `step_all_covariates_removed` - Numeric vector specifying the step at
#'  which all candidate covariates were deemed insignificant and removed. If the
#'  input data include competing risks, the vector records this step separately
#'  for each event. If the value is 0 for a given event, it means that at least
#'  one covariate survived the pre-selection through the final step of the
#'  procedure.
#'  * Per-event decision data frame(s) — one data frame per event (e.g., with
#'  competing risks one will get two: one for the event of interest and one for
#'  the competing event). The data frame name starts with `final_dec_pooled`
#'  (default) or `final_dec_per_variable`, depending on the method parameter.
#'  Contents depend on `to_show` and include: p-values and decision columns
#'  (names end with `_confirmed`) indicating whether a covariate is deemed
#'  informative at the final step of the procedure: 1 = informative, 0 = not
#'  informative. If all covariates were dropped during pre-selection (i.e., none
#'  reached the final step), all p-values are `NA` and all decisions are 0.
#'  * Per-event data frame named with the prefix `vimp_history`. If
#'  `save_vimp_history` is set to `"all"` or `"last"`, this data frame contains
#'  the VIMPs of covariates and their shadows from the final step of the
#'  procedure for the event indicated by the object’s name. If
#'  `save_vimp_history` is `"none"`, this element is `NULL`.
#'  * `time_elapsed` - list containing the run time of each step and the total
#'   time taken to execute the code.
#'  * `pre_selection` -  list in which the results of the pre-selection are
#'   stored. The exact form of this element depends on the chosen value of the
#'   `save_vimp_history` parameter.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @examples
#' data(veteran, package = "randomForestSRC")
#'
#' # When working with real data, use higher values for the niters and num.trees
#' # parameters --> here these parameters are set to small values to reduce the
#' # runtime.
#'
#' # Standard use - sruvival data without missing values or competing risks
#' out_veteran1 <- shadow_vimp_survival(
#'   data = veteran, time_column = "time", status_column = "status",
#'   niters = c(10, 20, 30), num.trees = 10
#' )
#'
#' # Save variable importance measures only from the final step of the
#' # procedure
#' out_veteran2 <- shadow_vimp_survival(
#'   data = veteran, time_column = "time", status_column = "status",
#'   niters = c(10, 20, 30), save_vimp_history = "last", num.trees = 10
#' )
#'
#' # Print unadjusted and FDR-adjusted p-values together with the corresponding
#' # decisions
#' out_veteran3 <- shadow_vimp_survival(
#'   data = veteran, time_column = "time", status_column = "status",
#'   niters = c(10, 20, 30), to_show = "FDR", num.trees = 10
#' )
#'
#' # Use per-variable p-values to decide in the final step whether a covariate
#' # is informative or not. Note that pooled p-values are always used in the
#' # pre-selection (first two steps).
#' out_veteran4 <- shadow_vimp_survival(
#'   data = veteran, time_column = "time", status_column = "status",
#'   niters = c(10, 20, 30), method = "per_variable", num.trees = 10
#' )
#'
#' # Perform variable selection in a single step, without a pre-selection phase
#' out_veteran5 <- shadow_vimp_survival(
#'   data = veteran, time_column = "time", status_column = "status",
#'   alphas = c(0.05), niters = c(30), num.trees = 10
#' )
#'
#' # shadowVIMP on survival data with missing values in the predictor
#' # For illustration, we create an artificial covariate with missing values for
#' # three subjects. By default shadow_vimp_survival(), excludes rows with
#' # missing values in any covariate.
#' new_var <-  sample.int(8, nrow(veteran) - 3, replace = TRUE)
#' veteran_na <- cbind(veteran, var1 = c(new_var, rep(NA, 3)))
#'
#' out_veteran_na <- shadow_vimp_survival(data = veteran_na,
#'  time_column = "time", status_column = "status",
#'  niters = c(5, 10, 12), num.trees = 10)
#'
#' # It is also possible to impute missing values
#' out_veteran_impute_na <- shadow_vimp_survival(data = veteran_na,
#'  time_column = "time", status_column = "status",
#'  na.action = "na.impute", niters = c(5, 10, 12), num.trees = 10)
#'
#' # shadowVIMP on survival data with competing risk
#' data(wihs, package = "randomForestSRC")
#' out_wihs <- shadow_vimp_survival(data = wihs,
#' time_column = "time", status_column = "status",
#' niters = c(5, 10, 12), num.trees = 10)
shadow_vimp_survival <- function(alphas = c(0.3, 0.10, 0.05),
                                 niters = c(30, 120, 1500),
                                 data,
                                 time_column,
                                 status_column,
                                 num.trees = max(2 * (ncol(data) - 1), 10000),
                                 importance = "permute",
                                 save_vimp_history = c("all", "last", "none"),
                                 to_show = c("FWER", "FDR", "unadjusted"),
                                 method = c("pooled", "per_variable"),
                                 na.action = c("na.omit", "na.impute"),
                                 ...) {
  cl <- match.call()
  cl[[1]] <- as.name("shadow_vimp_survival")

  # Check if there is the same number of alpha and niters parameters
  if (length(alphas) != length(niters)) {
    stop("`alphas` and `niters` must have the same length!")
  }
  if (is.unsorted(rev(alphas)) | any(alphas <= 0 | any(alphas >= 1))) {
    stop("Alphas must be in descending order. All alphas must be greater than 0 and less than 1.")
  }

  # Ensure save_vimp_history, to_show and method are one of the allowed values
  save_vimp_history <- match.arg(save_vimp_history)
  method <- match.arg(method)
  to_show <- match.arg(to_show)

  # Capture the name of the data object
  data_name <- deparse(substitute(data))

  # Capture the number of covariates (needed for correction of p-values)
  # ncol(data)-2 because we have time and status inidcator
  init_num_vars <- ncol(data) - 2

  replicate <- NULL
  vimp_sim_presel <- list()
  result_from_previous_step_bool <- list()

  for (j in 1:length(alphas)) {
    start_time <- Sys.time()
    message("alpha ", alphas[j], " \n")

    # Simulate VIMPs
    if (j == 1) {
      # run on full data
      vimpermsim <- vim_perm_sim_survival(
        data = data,
        time_column = time_column,
        status_column = status_column,
        niters = niters[j],
        importance = importance,
        num.trees = num.trees,
        data_name = data_name,
        na.action = na.action,
        ...
      )
      # Save event names
      event_names <- names(vimpermsim)
      # For the first iteration of algorithm, none of the results are taken from the previous step
      result_from_previous_step_bool <- setNames(as.list(rep(FALSE, length(event_names))), event_names)
    } else {
      # We are in at least 2nd step of pre-selection:
      for (event in event_names) {
        # For a given event, get the number of variables that survived previous pre-selection step
        no_vars_prev <- length(replicate[[j - 1]]$variables_remaining_for_replicate_pooled[[event]]) == 0
        result_from_previous_step_bool[[event]] <- no_vars_prev

        if (no_vars_prev) {
          next
        }

        # run only for events that still have some variables left after previous steps of pre-selection
        vimp_sim_presel[[event]] <- vim_perm_sim_survival(
          data = data %>%
            dplyr::select(dplyr::all_of(c(
              replicate[[j - 1]]$variables_remaining_for_replicate_pooled[[event]],
              status_column,
              time_column
            ))),
          time_column = time_column,
          status_column = status_column,
          niters = niters[j],
          importance = importance,
          num.trees = num.trees,
          data_name = data_name,
          na.action = na.action,
          ...
        )[[event]]
      }
    }

    # add test results only for events for which we simulated VIMPs in this step
    events_to_update <- names(result_from_previous_step_bool)[!unlist(result_from_previous_step_bool)]
    any_updates <- length(events_to_update) > 0

    if (any_updates) {
      if (j == 1) {
        # add_test_results to all events:
        vimpermsim_updated <- add_test_results_survival(
          vimpermsim,
          alpha = alphas[j],
          init_num_vars = init_num_vars,
          to_show = to_show
        )
      } else {
        # Apply add_test_results_survival only to the events for which we simulated VIMPs in the current step
        vimp_input <- vimp_sim_presel[events_to_update]
        vimpermsim_updated <- add_test_results_survival(
          vimp_input,
          alpha = alphas[j],
          init_num_vars = init_num_vars,
          to_show = to_show
        )
      }

      # Per event - make a list of remaining variables
      variables_remaining_for_replicate_pooled <- list()
      for (event in event_names) {
        if (isTRUE(result_from_previous_step_bool[[event]])) {
          # For the considered event no variables survived until now
          variables_remaining_for_replicate_pooled[[event]] <- replicate[[j - 1]]$variables_remaining_for_replicate_pooled[[event]]
          message(paste(event, ": no variables available at step", j, ".\nShowing previous-step results."))
        } else {
          # Get the variables that survived this step
          variables_remaining_for_replicate_pooled[[event]] <- vimpermsim_updated$test_res_pooled[[event]] %>%
            dplyr::filter(.data[["quantile_pooled"]] >= 1 - alphas[j]) %>%
            dplyr::select("varname") %>%
            unlist() %>%
            unname()
          message(event, "- variables remaining: ", length(variables_remaining_for_replicate_pooled[[event]]))
        }
      }

      # Start from previous (or first-step) object and update only the events we recomputed
      if (j == 1) {
        vimpermsim_to_store <- vimpermsim_updated
      } else {
        vimpermsim_to_store <- replicate[[j - 1]]$vimpermsim

        for (event in events_to_update) {
          # Update VIMPs and test results for the event for which we really computed this info in this step
          vimpermsim_to_store$test_res_pooled[[event]] <- vimpermsim_updated$test_res_pooled[[event]]
          vimpermsim_to_store$test_res_per_variable[[event]] <- vimpermsim_updated$test_res_per_variable[[event]]
          vimpermsim_to_store[[event]] <- vimpermsim_updated[[event]]
        }
      }

      end_time <- Sys.time()
      replicate[[j]] <- list(
        "alpha" = alphas[j],
        "result_taken_from_previous_step" = result_from_previous_step_bool,
        "vimpermsim" = vimpermsim_to_store,
        "variables_remaining_for_replicate_pooled" = variables_remaining_for_replicate_pooled,
        "time_elapsed" = difftime(end_time, start_time, units = "mins")
      )
    } else {
      # no updates for any event at this step: carry over the whole previous result
      end_time <- Sys.time()
      last_res <- replicate[[j - 1]]
      last_res$result_taken_from_previous_step <- result_from_previous_step_bool
      last_res$time_elapsed <- difftime(end_time, start_time, units = "mins")
      last_res$warning <- "No important variables remained for any event at this step. Results taken from previous step."
      replicate[[j]] <- last_res
    }
  }

  # time elapsed
  time <- list()
  for (i in 1:length(replicate)) {
    step_name <- paste0("step_", i)
    time[[step_name]] <- replicate[[i]]$time_elapsed %>%
      as.numeric()
  }

  time[["total_time_mins"]] <- time %>%
    unlist() %>%
    sum()

  # Clean vimp history if "none" or "last" option is selected
  for (event in event_names) {
    if (save_vimp_history == "none") {
      for (i in 1:length(replicate)) {
        replicate[[i]]$vimpermsim[[event]] <- NULL
      }
    } else if (save_vimp_history == "last") {
      for (i in 1:(length(replicate) - 1)) {
        replicate[[i]]$vimpermsim[[event]] <- NULL
      }
    }
  }

  # If the results were taken from the previous step, indicate which one
  vecs <- lapply(replicate, function(x) unlist(x$result_taken_from_previous_step, use.names = TRUE))
  from_prev_step <- do.call(cbind, vecs)

  # Named integer vector indicating for each event from which step the results were taken
  flag_from_prev_step <- apply(from_prev_step, 1, function(row) {
    idx <- which(row)
    if (length(idx) == 0) 0 else as.integer(idx[1] - 1)
  })

  # If pre-selection has been done - save its results
  if (length(alphas) > 1) {
    # Create pre_selection list storing results of the pre-selection steps
    pre_selection <- list()

    for (i in 1:(length(alphas) - 1)) {
      step_name <- paste0("step_", i)

      for (event in event_names) {
        name <- paste0("vimp_history_", event)
        pre_selection[[step_name]][[name]] <- replicate[[i]]$vimpermsim[[event]]

        name_dec <- paste0("decision_pooled_", event)
        if (flag_from_prev_step[[event]] == 0) {
          # Results were not taken from any previous step - some covariates survived until the end of the procedure
          pre_selection[[step_name]][[name_dec]] <- replicate[[i]]$vimpermsim$test_res_pooled[[event]]
        } else {
          pre_selection[[step_name]][[name_dec]] <- replicate[[i]]$vimpermsim$test_res_pooled[[event]] %>%
            mutate(
              across(starts_with("p_"), ~NA),
              across(starts_with("quantile_"), ~NA),
              across(ends_with("_confirmed"), ~0)
            )
        }
      }
      pre_selection[[step_name]][["alpha"]] <- alphas[i]
    }
  }

  # Create list storing results of the final step
  last_idx <- length(alphas)
  output <- list()

  for (event in event_names) {
    name <- paste0("vimp_history_", event)
    output[[name]] <- replicate[[last_idx]]$vimpermsim[[event]]
  }

  output[["alpha"]] <- alphas
  output[["step_all_covariates_removed"]] <- flag_from_prev_step
  output[["time"]] <- time
  output[["pre_selection"]] <- if (length(alphas) > 1) pre_selection else NULL
  output[["call"]] <- cl

  for (event in event_names) {
    if (flag_from_prev_step[[event]] == 0) {
      # If non-zero number of covariates survived until the last step of the procedure, report results from the last step
      if (method == "pooled") {
        final_dec <- replicate[[last_idx]]$vimpermsim$test_res_pooled[[event]] %>%
          select(-c("quantile_pooled"))
      } else {
        final_dec <- replicate[[last_idx]]$vimpermsim$test_res_per_variable[[event]] %>%
          select(-c("quantile_per_variable"))
      }
    } else {
      # None of the covariates survived until the last step of the procedure
      # In the final_dec output the user gets the list of all covariats, as p-values - NA,
      # In Type1_confirmed, etc. columns - only zeros
      final_dec <- replicate[[1]]$vimpermsim$test_res_pooled[[event]] %>%
        select(-c("quantile_pooled")) %>%
        mutate(
          across(starts_with("p_"), ~NA),
          across(ends_with("_confirmed"), ~0)
        )
    }

    sublist_name <- paste0("final_dec_", method, "_", event)
    output[[sublist_name]] <- final_dec
  }

  class(output) <- "shadow_vimp_survival"
  return(output)
}
