# Help Functions----------------------------------------------------------------

check_formula <- function(formula, data){
  checkmate::assert_formula(formula)
  
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
}


check_ransaclm_inputs <- function(formula, data, error_threshold,
                                  inlier_threshold, iterations,
                                  cores, seed) {
  checkmate::assert_data_frame(data)
  checkmate::assert_number(error_threshold)
  checkmate::assert_count(inlier_threshold)
  checkmate::assert_count(iterations)
  checkmate::assert_count(cores, null.ok = TRUE)
  checkmate::assert_number(seed, null.ok = TRUE)
  
  check_formula(formula, data)

}


stop_if_categorical <- function(model_matrix) {
  numb_col <- ncol(model_matrix)
  for (i in seq_len(numb_col)) {
    unique_values <- unique(model_matrix[, i])
    if (length(unique_values) == 2) {
      stop("Sorry, but this version does not work for categorical variables.
           Please remove them.")
    }
  }
}


data_dimension <- function(formula, data) {
  dummy_lm <- lm(formula, data = data)

  model_matrix <- model.matrix(dummy_lm)

  min_subset_size <- ncol(model_matrix)
  numb_observations <- nrow(model_matrix)

  if (numb_observations <= 2) {
    stop("Not enough observations. With one or two observations it is not really
         meaningful to talk about outliers.")
  }
  
  if (min_subset_size - 1 > numb_observations) {
    stop("Not enough observations. More observations than variables are needed!")
  }
  
  stop_if_categorical(model_matrix)

  c("min_size" = min_subset_size, "Observations" = numb_observations)
}


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

    mean_squared_error <- mean(new_model$residuals^2)
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


best_model <- function(error_inlier_list, data, numb_observations = nrow(data)) {
  which_min_error <- which_min_error(error_inlier_list)

  model <- error_inlier_list[[which_min_error]]$model

  inlier_subset <- error_inlier_list[[which_min_error]]$inlier

  data$.consensus_set <- is.element(seq_len(numb_observations), inlier_subset)

  list("data" = data, "model" = model)
}


# Wurde nicht verwendet, da bei inlier_fraction = 0.3, n_obs = 100 und
# successs_prob zB 0.8, das Ergebnis ca 3*10^58 wäre und für soviel
# Iterationen hat niemand Zeit.
numb_of_iterations <- function(n_obs, inlier_fraction, success_prob) {
  as.integer(
    log(1 - success_prob) / log(1 - inlier_fraction^n_obs)
  )
}


# ransaclm----------------------------------------------------------------------
ransaclm <- function(formula, data, error_threshold, inlier_threshold,
                     iterations = 500,
                     cores = NULL, seed = NULL) {
  check_ransaclm_inputs(
    formula, data, error_threshold,
    inlier_threshold, iterations, cores, seed
  )

  if (sum(is.na(data)) > 0) {
    data <- na.omit(data)
    warning("Removed at least one observation with missing data!")
  }

  dimension <- data_dimension(formula, data)

  min_subset_size <- dimension["min_size"]
  numb_observations <- dimension["Observations"]


  observed_response <- data[, as.character(formula)[2]]

  set.seed(seed = seed)

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
  best_model(error_inlier_list,
    data = data, numb_observations = numb_observations
  )
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


# Own tests---------------------------------------------------------------------


testthat::test_that("Warning for Missing data", {
  data_na <- data_simple
  data_na[1, 1] <- NA


  testthat::expect_warning(
    ransaclm(y ~ . - inlier,
      data = data_na, error_threshold = 2,
      inlier_threshold = 50, seed = 20171111
    )
  )
})


testthat::test_that("Error for categorical data", {
  testthat::expect_error(
    ransaclm(y ~ .,
      data = data_simple, error_threshold = 2,
      inlier_threshold = 50, seed = 20171111
    )
  )

  data_cat <- data_simple
  data_cat$abc <- sample(c("a", "b", "c"), size = 100, replace = TRUE)

  testthat::expect_error(
    ransaclm(y ~ x + abc,
      data = data_simple, error_threshold = 2,
      inlier_threshold = 50, seed = 20171111
    )
  )
})


testthat::test_that("Error for nonsense formula and not enough data", {
  testthat::expect_error(
    ransaclm(y~gibtesnicht,
             data = data_simple, error_threshold = 2,
             inlier_threshold = 50, seed = 20171111
    )
  )
  
  data_less <- data_simple[1:2,]
  data_p_bigger_n <- data_simple[1:3,]
  data_p_bigger_n$a <- rnorm(3)
  data_p_bigger_n$b <- rnorm(3)
  data_p_bigger_n$c <- rnorm(3)
  
  testthat::expect_error(
    ransaclm(y ~ x ,
             data = data_less, error_threshold = 2,
             inlier_threshold = 50, seed = 20171111
    )
  )
  testthat::expect_error(
    ransaclm(y ~ x + a + b + c,
             data = data_p_bigger_n, error_threshold = 2,
             inlier_threshold = 50, seed = 20171111
    )
  )
  
})


testthat::test_that("Works normally", {
  # No intercept
  ransaclm(y ~ . - inlier - 1,
    data = data_simple, error_threshold = 2,
    inlier_threshold = 50, seed = 20171111
  )
  # Other formula notation
  ransaclm(y ~ x,
    data = data_simple, error_threshold = 2,
    inlier_threshold = 50, seed = 20171111
  )
  
  # Add dimension
  data2 <- data_simple
  data2$z <- rnorm(100)
  
  ransaclm(y ~ x + z,
           data = data2, error_threshold = 2,
           inlier_threshold = 50, seed = 20171111
  )
  
  ransaclm(y ~ . - inlier,
           data = data2, error_threshold = 2,
           inlier_threshold = 50, seed = 20171111
  )
  
})