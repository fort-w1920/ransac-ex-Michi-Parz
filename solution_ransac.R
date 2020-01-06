# Help Functions----------------------------------------------------------------

subset_error <- function(subsample, data, formula, observed_response,
                         error_threshold, inlier_threshold) {
  new_model <- NULL

  mean_squared_error <- Inf

  subset <- data[subsample, ]

  linear_model <- lm(formula = formula, data = subset)

  prediction <- predict.lm(linear_model, data)

  absolute_distance <- abs(prediction - observed_response)

  inliers_subset <- which(absolute_distance <= error_threshold)



  if (length(inliers_subset) >= inlier_threshold) {
    inlier_data <- data[inliers_subset, ]

    new_model <- lm(formula = formula, data = inlier_data)

    mean_squared_error <- mean(new_model$residuals ^ 2)
  }


  list(
    "error" = mean_squared_error,
    "inlier" = inliers_subset,
    "model" = new_model
  )
}


which_min_error <- function(error_inlier_list) {
  unlist_error_inlier <- unlist(error_inlier_list)

  which_error <- which(names(unlist_error_inlier) == "error")

  unlist_error <- unlist_error_inlier[which_error]

  which.min(unlist_error)
}


numb_of_iterations <- function(n_obs, inlier_fraction, success_prob) {
  as.integer(
    log(1 - success_prob) / log(1 - inlier_fraction ^ n_obs)
  )
}


check_ransaclm_inputs <- function(formula, data, error_threshold,
                                  inlier_threshold, iterations,
                                  cores, seed) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_number(error_threshold)
  checkmate::assert_count(inlier_threshold)
  checkmate::assert_count(iterations)
  checkmate::assert_count(cores, null.ok = TRUE)
  checkmate::assert_number(seed, null.ok = TRUE)


  formula_variables <- all.vars(formula)

  formula_dot_ex <- formula_variables[2] == "."


  col_names <- names(data)

  variables_exists <- is.element(
    formula_variables,
    col_names
  )

  if (!all(variables_exists) && !formula_dot_ex) {
    stop("At least one variable, used in formula, is not in your data!")
  }

  dummy_lm <- lm(formula, data = data)

  #data_classes <- attr(summary(dummy_lm)$terms, "dataClasses")


  dummy_model <- dummy_lm$model

  if (nrow(dummy_model) < ncol(dummy_model)) {
    stop("More observations than covariables are needed!")
  }
}

# ransaclm----------------------------------------------------------------------
ransaclm <- function(formula, data, error_threshold, inlier_threshold,
                     iterations = floor(nrow(data) / 3),
                     cores = NULL, seed = NULL) {
  check_ransaclm_inputs(
    formula, data, error_threshold,
    inlier_threshold, iterations, cores, seed
  )

  set.seed(seed = seed)

  min_subset_size <- ncol(data)
  numb_observations <- nrow(data)
  observed_response <- data[, as.character(formula)[2]]

  subsamples <- replicate(
    iterations,
    sample(seq_len(numb_observations),
      size = min_subset_size
    ),
    simplify = FALSE
  )


  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1L
  }

  cluster <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl = cluster))


  error_inlier_list <- parallel::parLapply(
    cl = cluster,
    X = subsamples,
    fun = subset_error,
    data = data,
    formula = formula,
    observed_response = observed_response,
    error_threshold = error_threshold,
    inlier_threshold = inlier_threshold
  )


  # Return the best Model
  which_min_error <- which_min_error(error_inlier_list)

  model <- error_inlier_list[[which_min_error]]$model

  inlier_subset <- error_inlier_list[[which_min_error]]$inlier

  data$.consensus_set <- is.element(seq_len(numb_observations), inlier_subset)

  list("data" = data, "model" = model)
}


# load ransac-utils.R-----------------------------------------------------------
source("ransac-utils.R")


# ransac-example----------------------------------------------------------------
set.seed(1874374111)
data_simple <- make_ransac_data(n_obs = 100, n_coef = 1, inlier_fraction = 0.7)

# univariate example:
ransac_simple <- ransaclm(y ~ . - inlier,
  data = data_simple, error_threshold = 2,
  inlier_threshold = 50, seed = 20171111
)
validate_ransac(ransac_simple)
