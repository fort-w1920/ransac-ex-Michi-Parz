
ransaclm <- function(formula, data, error_threshold, inlier_threshold,
                     iterations) {

  min_subset_size <- ncol(data)
  numb_observations <- nrow(data)
  observed_response <- data[,as.character(formula)[2]]
  
  lowest_error <- Inf

  for (i in seq_len(iterations)) { # Kann man auch paralellisieren aber erstmal so
    subsamples <- sample(seq_len(numb_observations), size = min_subset_size)
    subset <- data[subsamples,]
    
    linear_model <- lm(formula = formula, data = subset)
    
    prediction <- predict.lm(linear_model, data)
    
    absolute_distance <- abs(prediction - observed_response) 
    
    inliers_subset <- which(absolute_distance <= error_threshold)
    
    if (inliers_subset >= inlier_threshold) {
      inlier_data <- data[inliers_subset,]
      
      new_model <- lm(formula = formula, data = inlier_data)
      
      mean_squared_error <- mean(new_model$residuals^2)
      
      if (mean_squared_error < lowest_error) {
        lowest_error <- mean_squared_error
        
        result <- list(new_model, inlier_data)
        
      }
      
      
    }
    
  }



  result
}


#Rumprobieren
x1 <- rnorm(100)
x2 <- rnorm(100)
eps <- rnorm(100, sd = 0.2)

y <- x1 + x2 + eps

df <- data.frame(y,x1,x2)


formular <- y~.
lm <- lm(formular,data = df)

summary(lm)

lm$effects

summary_lm <- summary(lm)

summary_lm
lm

?lm

