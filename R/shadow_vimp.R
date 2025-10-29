#' Select influential covariates by using random forests with multiple testing
#' control
#'
#' `shadow_vimp()` performs variable selection and determines whether each
#' covariate is influential based on unadjusted, FDR-adjusted, and FWER-adjusted
#' p-values.
#'
#'
#' The `shadow_vimp()` function by default performs variable selection in
#' multiple steps. Initially, it prunes the set of predictors using a relaxed
#' (higher) alpha threshold in a pre-selection stage. Variables that pass this
#' stage then undergo a final evaluation using the target (lower) alpha
#' threshold and more iterations. This stepwise approach distinguishes
#' informative from uninformative covariates based on their VIMPs and enhances
#' computational efficiency. The user can also perform variable selection in a
#' single step, without a pre-selection phase.
#'
#' @param alphas Numeric vector, significance level values for each step of the
#'   procedure, default `c(0.3, 0.10, 0.05)`.
#' @param niters  Numeric vector, number of permutations to be performed in each
#'   step of the procedure, default `c(30, 120, 1500)`.
#' @param data Input data frame.
#' @param outcome_var Character, name of the column containing the outcome
#'   variable.
#' @param num.threads Numeric. The number of threads used by [ranger::ranger()]
#'   for parallel tree building. If `NULL` (the default), half of the available
#'   CPU threads are used (this is the default behaviour in `shadow_vimp()`,
#'   which is different from the default in [ranger::ranger()]). See the
#'   [ranger::ranger()] documentation for more details.
#' @param to_show Character, one of `"FWER"`, `"FDR"` or `"unadjusted"`.
#'  * `"FWER"` (the default) - the output includes unadjusted,
#'   Benjamini-Hochberg (FDR) and Holm (FWER) adjusted p-values together with
#'   the decision whether the variable is significant or not (1 - significant, 0
#'   means not significant) according to the chosen criterium.
#'  * `"FDR"` - the output includes both unadjusted and FDR adjusted p-values
#'  along with the decision.
#'  * `"unadjusted:` - the output contains only raw, unadjusted p-values
#'   together with the decision.
#' @param num.trees Numeric, number of trees. Passed to [ranger::ranger()],
#'   default is `max(2 * (ncol(data) - 1), 10000)`.
#' @param importance Character, the type of variable importance to be calculated
#'   for each variable. Argument passed to [ranger::ranger()], default is
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
#' @param ... Additional parameters passed to [ranger::ranger()].
#' @return Object of the class "shadow_vimp" with the following entries:
#'  * `call` - the call formula used to generate the output.
#'  * `alpha` - numeric, significance level used in the algorithm.
#'  * `step_all_covariates_removed` - integer. If > 0, the step number at which
#'   all candidate covariates were deemed insignificant and removed. If 0, at
#'   least one covariate survived the pre-selection until the last step of the
#'   procedure.
#'  * `final_dec_pooled` (the default) or `final_dec_per_variable` -  a data
#'   frame that contains, depending on the specified value of the `to_show`
#'   parameter, p-values and corresponding decisions (in columns with names
#'   ending in `confirmed`) if the variable is deemed informative at the final
#'   step of the procedure: 1 = covariate considered informative in the last
#'   step; 0 = not informative. If all covariates were dropped in the
#'   pre-selection, i.e. none reached the final step, then all p-values are NA
#'   and all decisions are set to 0.
#'  * `vimp_history`- if `save_vimp_history` is set to `"all"` or `"last"` then
#'   it is a data frame with VIMPs of covariates and their shadows from the last
#'   step of the procedure. If `save_vimp_history` is set to `"none"`, then it
#'   is `NULL`.
#'  * `time_elapsed` - list containing the runtime of each step and the total
#'   time taken to execute the code.
#'  * `pre_selection` -  list in which the results of the pre-selection are
#'   stored. The exact form of this element depends on the chosen value of the
#'   `save_vimp_history` parameter.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @examples
#' data(mtcars)
#'
#' # When working with real data, use higher values for the niters and num.trees
#' # parameters --> here these parameters are set to small values to reduce the
#' # runtime.
#'
#' # Function to make sure proper number of cores is specified
#' safe_num_threads <- function(n) {
#'   available <- parallel::detectCores()
#'   if (n > available) available else n
#' }
#'
#' # Standard use
#' out1 <- shadow_vimp(
#'   data = mtcars, outcome_var = "vs",
#'   niters = c(10, 20, 30), num.trees = 30, num.threads = safe_num_threads(1)
#' )
#'
#' \donttest{
#' # `num.threads` sets the number of threads for multithreading in
#' # `ranger::ranger`. By default, the `shadow_vimp` function uses half the
#' # available CPU threads.
#' out2 <- shadow_vimp(
#'   data = mtcars, outcome_var = "vs",
#'   niters = c(10, 20, 30), num.threads = safe_num_threads(2),
#'   num.trees = 30
#' )
#'
#' # Save variable importance measures only from the final step of the
#' # procedure
#' out4 <- shadow_vimp(
#'   data = mtcars, outcome_var = "vs",
#'   niters = c(10, 20, 30), save_vimp_history = "last", num.trees = 30,
#'   num.threads = safe_num_threads(1)
#' )
#'
#' # Print unadjusted and FDR-adjusted p-values together with the corresponding
#' # decisions
#' out5 <- shadow_vimp(
#'   data = mtcars, outcome_var = "vs",
#'   niters = c(10, 20, 30), to_show = "FDR", num.trees = 30,
#'   num.threads = safe_num_threads(1)
#' )
#'
#' # Use per-variable p-values to decide in the final step whether a covariate
#' # is informative or not. Note that pooled p-values are always used in the
#' # pre-selection (first two steps).
#' out6 <- shadow_vimp(
#'   data = mtcars, outcome_var = "vs",
#'   niters = c(10, 20, 30), method = "per_variable", num.trees = 30,
#'   num.threads = safe_num_threads(1)
#' )
#'
#' # Perform variable selection in a single step, without a pre-selection phase
#' out7 <- shadow_vimp(
#'   data = mtcars, outcome_var = "vs", alphas = c(0.05),
#'   niters = c(30), num.trees = 30,
#'   num.threads = safe_num_threads(1)
#' )
#' }
shadow_vimp <- function(alphas = c(0.3, 0.10, 0.05),
                        niters = c(30, 120, 1500),
                        data,
                        outcome_var, # y,
                        num.trees = max(2 * (ncol(data) - 1), 10000),
                        num.threads = NULL,
                        importance = "permutation",
                        save_vimp_history = c("all", "last", "none"),
                        to_show = c("FWER", "FDR", "unadjusted"),
                        method = c("pooled", "per_variable"),
                        ...) {
  cl <- match.call()
  cl[[1]] <- as.name("shadow_vimp")

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
  init_num_vars <- ncol(data) - 1

  # for loop over alphas
  replicate <- NULL

  for (j in 1:length(alphas)) {
    # runtime start
    start_time <- Sys.time()

    message("alpha ", alphas[j], " \n")

    if (j > 1 && length(replicate[[j - 1]]$variables_remaining_for_replicate_pooled) == 0) {
      result_from_previous_step_bool <- TRUE
      warning("None of the variables have passed the pre-selection process to the final step.")
    } else {
      # Warnings concerning the number of iterations if the user specified too small number
      if (j < length(alphas)) {
        # We are in the pre-selection phase, where always pooled method is used
        if (j > 1) {
          # Number of available covariates = number of covariates that survived pre-selection in the previous step
          p <- data %>%
            select(all_of(c(replicate[[j - 1]]$variables_remaining_for_replicate_pooled))) %>%
            ncol()
        } else {
          # We are in the first step of the procedure so number of available variables = total number of variables
          p <- init_num_vars
        }

        # Theoretically the smallest p-value we can get:
        p_val_min <- 1 / (niters[j] * p)

        if (p_val_min > alphas[j] / init_num_vars) {
          warning("Not enough iterations for any positives after FDR/FWER adjustment.\n Increase the number of iterations in the pre-selection phase to get reliable results.")
        }
      } else {
        # We are in the final step of the procedure - decision can be made using pooled or per_variable approach
        # Number of available variables in the last step:
        p <- data %>%
          select(all_of(c(replicate[[j - 1]]$variables_remaining_for_replicate_pooled))) %>%
          ncol()

        if (method == "pooled") {
          # Theoretically the smallest p-value we can get when using pooled approach:
          p_val_min <- 1 / (niters[j] * p)

          if (p_val_min > alphas[j] / init_num_vars) {
            warning("Not enough iterations for any positives after FDR/FWER adjustment.\n Increase the number of iterations in the final step to get reliable results.")
          }
        } else {
          # Theoretically the smallest p-value we can get when using per_variable approach:
          p_val_min <- 1 / (niters[j])

          if (p_val_min > alphas[j] / init_num_vars) {
            warning("Not enough iterations for any positives after FDR/FWER adjustment.\n Increase the number of iterations in the final step to get reliable results.")
          }
        }
      }

      # run algorithm
      vimpermsim <- vim_perm_sim(
        data = if (j > 1) {
          data %>%
            select(all_of(c(replicate[[j - 1]]$variables_remaining_for_replicate_pooled, outcome_var)))
        } else {
          data
        }, # pooled pre selection
        outcome_var = outcome_var, # y
        niters = niters[j],
        importance = importance,
        num.threads = num.threads,
        num.trees = num.trees,
        data_name = data_name,
        ...
      )

      result_from_previous_step_bool <- FALSE
    }

    # If there are some variables remaining after the pre-selection
    if (result_from_previous_step_bool == FALSE) {
      vimpermsim <- add_test_results(vimpermsim,
        alpha = alphas[j],
        init_num_vars = init_num_vars,
        to_show = to_show
      )

      variables_remaining_for_replicate_pooled <- vimpermsim$test_results$pooled %>%
        filter(.data[["quantile_pooled"]] >= 1 - alphas[j]) %>%
        select("varname") %>%
        unlist() %>%
        unname()

      message("Variables remaining: ", length(variables_remaining_for_replicate_pooled), "\n")

      # runtime end
      end_time <- Sys.time()

      replicate[[j]] <- list(
        "alpha" = alphas[j],
        "result_taken_from_previous_step" = result_from_previous_step_bool,
        "vimpermsim" = vimpermsim,
        "variables_remaining_for_replicate_pooled" = variables_remaining_for_replicate_pooled,
        "time_elapsed" = difftime(end_time, start_time, units = "mins")
      )
    } else {
      # no variables survived the pre-selection process --> take results from previous step
      end_time <- Sys.time()
      last_res <- replicate[[j - 1]]
      last_res$result_taken_from_previous_step <- result_from_previous_step_bool
      last_res$time_elapsed <- difftime(end_time, start_time, units = "mins")
      last_res$warning <- "No important variables remained until the last step, results taken from previous step."
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
  if (save_vimp_history == "none") {
    for (i in 1:length(replicate)) {
      replicate[[i]]$vimpermsim$vim_simulated <- NULL
    }
  } else if (save_vimp_history == "last") {
    for (i in 1:(length(replicate) - 1)) {
      replicate[[i]]$vimpermsim$vim_simulated <- NULL
    }
  }

  # If the results were taken from the previous step, indicate which one
  idx_from_prev_step <- which(sapply(replicate, function(x) x$result_taken_from_previous_step))[1]
  flag_from_prev_step <- if_else(is.na(idx_from_prev_step), 0, idx_from_prev_step - 1)

  # If pre-selection has been done - save its results
  if (length(alphas) > 1) {
    # Create pre_selection list storing results of the pre-selection steps
    pre_selection <- list()

    for (i in 1:(length(alphas) - 1)) {
      step_name <- paste0("step_", i)
      pre_selection[[step_name]][["vimp_history"]] <- replicate[[i]]$vimpermsim$vim_simulated

      if (flag_from_prev_step == 0) {
        # Results were not taken from any previous step - some covariates survived until the end of the procedure
        pre_selection[[step_name]][["decision_pooled"]] <- replicate[[i]]$vimpermsim$test_results$pooled
      } else {
        pre_selection[[step_name]][["decision_pooled"]] <- replicate[[i]]$vimpermsim$test_results$pooled %>%
          mutate(
            across(starts_with("p_"), ~NA),
            across(starts_with("quantile_"), ~NA),
            across(ends_with("_confirmed"), ~0)
          )
      }

      pre_selection[[step_name]][["alpha"]] <- alphas[i]
      # pre_selection[[step_name]][["alpha"]] <- replicate[[i]]$alpha
      # pre_selection[[step_name]][["result_taken_from_previous_step"]] <- replicate[[i]]$vimpermsim$result_taken_from_previous_step
    }
  }

  # Create list storing results of the final step
  last_idx <- length(alphas)

  output <- list(
    "vimp_history" = replicate[[last_idx]]$vimpermsim$vim_simulated,
    "alpha" = alphas,
    "step_all_covariates_removed" = flag_from_prev_step,
    "time_elapsed" = time,
    "pre_selection" = if (length(alphas) > 1) pre_selection else NULL,
    "call" = cl
  )

  if (flag_from_prev_step == 0) {
    # If non-zero number of covariates survived until the last step of the procedure, report results from the last step
    if (method == "pooled") {
      final_dec <- replicate[[last_idx]]$vimpermsim$test_results$pooled %>% select(-.data[["quantile_pooled"]])
    } else {
      # if(flag_from_prev_step == 0){
      final_dec <- replicate[[last_idx]]$vimpermsim$test_results$per_variable %>% select(-.data[["quantile_per_variable"]])
      # }
    }
  } else {
    # None of the covariates survived until the last step of the procedure
    # In the final_dec output the user gets the list of all covariats, as p-values - NA,
    # In Type1_confirmed, etc. columns - only zeros
    final_dec <- replicate[[1]]$vimpermsim$test_results$pooled %>%
      select(-.data[["quantile_pooled"]]) %>%
      mutate(
        across(starts_with("p_"), ~NA),
        across(ends_with("_confirmed"), ~0)
      )
  }


  sublist_name <- paste0("final_dec_", method)
  output[[sublist_name]] <- final_dec

  class(output) <- "shadow_vimp"
  return(output)
}
