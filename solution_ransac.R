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

    mean_squared_error <- mean(new_model$residuals^2)
  }


  list(
    "error" = mean_squared_error,
    "inlier" = inliers_subset,
    "model" = new_model
  )
}


which_min_error <- function(error_inlier_list){
  unlist_error_inlier <- unlist(error_inlier_list)
  
  which_error <- which(names(unlist_error_inlier) == "error")
  
  unlist_error <- unlist_error_inlier[which_error]
  
  which.min(unlist_error)
}



check_ransaclm_inputs <- function(formula, data, error_threshold,
                                  inlier_threshold, iterations,
                                  cores, seed) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_number(error_threshold)
  checkmate::assert_count(inlier_threshold)
  checkmate::assert_number(iterations, lower = 1)
  checkmate::assert_count(cores, null.ok = TRUE)
  checkmate::assert_number(seed, null.ok = TRUE)
}

# ransaclm----------------------------------------------------------------------
ransaclm <- function(formula, data, error_threshold, inlier_threshold,
                     iterations = nrow(data) / 3, cores = NULL, seed = NULL) {
  check_ransaclm_inputs(
    formula, data, error_threshold,
    inlier_threshold, iterations, cores, seed
  )
  
  set.seed(seed = seed)

  min_subset_size <- ncol(data)
  numb_observations <- nrow(data)
  observed_response <- data[, as.character(formula)[2]]

  seq_iterations <- seq_len(iterations)

  # Subsamples
  subsamples <- lapply(seq_iterations, sample,
    x = seq_len(numb_observations),
    size = min_subset_size
  )
  
  # Make cluster
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1L
  }

  cluster <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl = cluster))

  # Error list
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


# Rumprobieren
x1 <- rnorm(100)
x2 <- rnorm(100)
eps <- rnorm(100, sd = 0.2)

y <- x1 + x2 + eps

df <- data.frame(y, x1, x2)


formular <- y ~ .
lm <- lm(formular, data = df)

summary(lm)

lm$effects

summary_lm <- summary(lm)

summary_lm
lm

?lm
