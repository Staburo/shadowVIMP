#' Compute permutation variable importance for predictors and their row-wise
#' shadows in the survival-data setting.
#'
#' `vim_perm_sim_survival()` calculates repeatedly (`niters` times) the variable
#' importance of the original values of the predictors and their row-wise
#' permuted shadows. Each shadow's variable importance is computed based on a
#' new permutation of the initial predictor values. This is the survival
#' equivalent of `vim_perm_sim()`.
#'
#' @param data Input data frame.
#' @param time_column Character, name of the column containing the event time.
#' @param status_column Character, name of the column containing the status
#'  indicator. The current implementation of the function supports the
#'  competing-risks setting.
#' @param niters Numeric, number of permutations of the initial predictor
#'   values, default is 100.
#' @param importance Character, the type of variable importance to be calculated
#'   for each independent variable. Argument passed to [randomForestSRC::rfsrc()],
#'   default is `permutation`.
#' @param num.trees Numeric, number of trees. Passed to [randomForestSRC::rfsrc()],
#'   default is `max(2 * (ncol(data) - 1), 10000)`.
#' @param data_name Character, name of the object passed as `data`. In
#'   `shadow_vimp()` it is set automatically.
#' @param na.action Character, one of `"na.omit"` or `"na.impute"`. Action to
#' take when the data contain `NA`. See [randomForestSRC::rfsrc()] and the
#' vignette on this package’s survival-data extension for details.
#' @param ... Additional parameters passed to [randomForestSRC::rfsrc()].
#' @return List containing `niters` variable importance values for both the
#'   original and row-wise permuted predictors.
#' @noRd
#' @import rlang dplyr
#' @importFrom magrittr %>%
#' @importFrom randomForestSRC rfsrc
#' @importFrom purrr map
#' @examples
#' # Standard survival data: Veterans' Administration Lung Cancer study data
#' data(veteran, package = "randomForestSRC")
#'
#' # When working with real data, increase num.trees value or keep the default.
#' # Here this parameter is set to a small value in order to reduce the run time.
#' # Standard use:
#' vps_out <- vim_perm_sim_survival(data = veteran,
#' time_column = "time",
#' status_column = "status",
#' niters = 30,
#' num.trees = 30)
vim_perm_sim_survival <- function(data,
                                  time_column,
                                  status_column,
                                  niters = 100,
                                  importance = "permute",
                                  num.trees = max(2 * (ncol(data) - 1), 10000),
                                  data_name = NULL,
                                  na.action = c("na.omit", "na.impute"),
                                  ...) {
  # Check if niters parameter has a correct format
  if (is.numeric(niters) == F) {
    stop("`niters` parameter must be numeric.")
  }

  # Capture the name of the object passed as data - needed for progress tracking
  if (is.null(data_name)) {
    data_name <- deparse(substitute(data))
  }

  # Check if passed data_name argument has a correct format
  if (length(data_name) > 1) {
    stop("`data_name` must be a character of length 1.")
  }

  if (!is.character(time_column) || !is.character(status_column) || length(time_column) != 1 || length(status_column) != 1) {
    stop("`time_column` and  `status_column` must be characters of length 1.")
  }

  p <- ncol(data) - 2
  n <- nrow(data)

  # Check the number of events & decide on the splitting rule
  num_events <- data |>
    pull(status_column) |>
    unique() |>
    length()
  if (num_events == 2) {
    # We are in the standard survival analysis setting
    split_rule <- "logrank"
  } else if (num_events > 2) {
    # We are in the competing risks setting
    split_rule <- "logrankCR"
  }

  # Splitting predictors
  predictors <- data %>% select(-all_of(c(status_column, time_column)))

  if (sum(grepl("_permuted$", names(predictors))) > 0) {
    stop("One or more variables ending with _permuted. Please rename them.")
  }

  # Creating permuted predictors
  predictors_p <- predictors[sample(1:n), , drop = FALSE]

  names(predictors_p) <- paste0(names(predictors_p), "_permuted")

  # Concatenating predictors, permuted predictors and label
  dt <- cbind.data.frame(
    predictors,
    predictors_p,
    data %>% select(all_of(c(status_column, time_column)))
  )

  # Prepare formula for RF
  # turn characters into a symbols
  time_sym <- ensym(time_column)
  status_sym <- ensym(status_column)

  formula_rf <- as.formula(
    expr(Surv(!!time_sym, !!status_sym) ~ .)
  )

  # rfsrc() has no analogue of respect.unordered.factors = "order" from ranger
  # Therefore I do the ordering of factors by hand in the same way as in ranger
  # --> for survival data factors are ordered by the median survival time
  to_convert <- dt %>%
    select(-all_of(c(status_column, time_column))) %>%
    select(where(is.factor)) %>%
    names()

  for (var in to_convert) {
    var_sym <- sym(var)

    order <- dt %>%
      group_by(!!var_sym) %>%
      summarise(
        median_time = median(!!time_sym, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(median_time) %>%
      pull(!!var_sym)

    dt <- dt %>%
      mutate(!!var_sym := factor(
        !!var_sym,
        levels  = order,
        ordered = TRUE
      ))
  }

  # Simulation
  vimp_sim <- NULL

  for (i in 1:niters) {
    if ((i %% 50 == 0) | (i == 1)) {
      message(paste0(
        format(Sys.time()), ": dataframe = ", data_name, " niters = ",
        niters, " num.trees = ", num.trees, ". Running step ", i, "\n"
      ))
    }

    # reshuffle row wise
    dt[, (ncol(predictors_p) + 1):(2 * ncol(predictors_p))] <- dt[sample(1:n), (ncol(predictors_p) + 1):(2 * ncol(predictors_p))]

    imp <- rfsrc(
      formula = formula_rf,
      data = dt,
      ntree = num.trees,
      splitrule = split_rule,
      importance = importance,
      na.action = na.action, # How to handle NAs
      samptype = "swr", # Sample with replacement
      save.memory = TRUE, # Set to save memory and speed up computation
      ...
    )$importance

    imp_df <- if (is.matrix(imp)) {
      as.data.frame(imp)
    } else {
      data.frame(event = imp, check.names = FALSE)
    }

    vimp_sim[[i]] <- imp_df
  }

  # Get the competing events names
  event_names <- colnames(vimp_sim[[1]])

  # Collect VIMPs separately for each of competing events
  results_sim <- list()
  for (event in event_names) {
    sim_event <- map(vimp_sim, ~ .x %>%
      select(all_of(event)) %>%
      t() %>%
      data.frame(row.names = NULL)) %>%
      bind_rows()
    results_sim[[event]] <- sim_event
  }

  # Check if for any event any covariate has always VIMP = 0
  for (event in event_names) {
    sd_shadow <- results_sim[[event]] %>%
      select(ends_with("_permuted")) %>%
      summarise(across(everything(), ~ sd(.x))) %>%
      pivot_longer(cols = everything(), names_to = "variable", values_to = "sd") %>%
      filter(sd == 0)

    if (nrow(sd_shadow) > 0) {
      warning("One or more shadow variables always have VIMP equal to zero.")
    }
  }

  return(results_sim)
}
