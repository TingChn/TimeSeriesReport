library(readr)
library(glmnet)
library(HDeconometrics)
library(ggplot2)
library(hdm)
library(dplyr)
library(lubridate)
library(tidyr)

#-----------------------------------------------------------------------------------------------------------
#Data processing 
#-----------------------------------------------------------------------------------------------------------
# Read the CSV file
raw <- read_delim("bk_data.csv", delim = ";", show_col_types = FALSE)

# Convert to a timetable-like structure
raw_tibble <- raw %>%
  mutate(Time = as.Date(raw[[1]], format = "%m.%d.%Y")) %>%
  select(-1) %>%
  relocate(Time)

# Extract metadata
t_code <- raw_tibble[1, -1] %>% as.character()
varlist <- colnames(raw_tibble)[-1]

# Remove the first row (metadata) from the data
data_tibble <- raw_tibble[-1, ]

# Convert numeric columns to appropriate data types
data_tibble <- data_tibble %>%
  mutate(across(-Time, ~ as.numeric(as.character(.))))

# Construct stationary series
#==========================================================================
data_trans <- data_tibble

for (ii in seq_along(varlist)) {
  col_name <- varlist[ii]
  transform_type <- t_code[ii] 
  
  if (transform_type == "fl") {
    if (any(data_trans[[col_name]] <= 0, na.rm = TRUE)) {
      # Handle non-positive values
      data_trans[[col_name]] <- c(NA, 100 * (data_trans[[col_name]][-1] / 
                                               data_trans[[col_name]][-nrow(data_trans)] - 1))
    } else {
      # Log-difference
      data_trans[[col_name]] <- c(NA, 100 * diff(log(data_trans[[col_name]])))
    }
  } else if (transform_type == "fd") {
    # First difference
    data_trans[[col_name]] <- c(NA, diff(data_trans[[col_name]]))
  }
  # If 'lv', do nothing (keep level)
}

# Choose final sample and remove series with missings
#==========================================================================
est_start <- as.Date("1965-01-01")
est_end <- as.Date("2019-12-01")

# Filter rows by time range
data_trans <- data_trans %>%
  filter(Time >= est_start & Time <= est_end)

# Remove columns with missing values
data_in <- data_trans %>%
  select(where(~ all(!is.na(.x))))

# Convert to matrix for further analysis
data <- as.matrix(data_in %>% select(-Time))

# Dimensions of data
T <- nrow(data)
N <- ncol(data)

# Standardize the input data
#==========================================================================
data_mean <- colMeans(data, na.rm = TRUE)
data_sd <- apply(data, 2, sd, na.rm = TRUE)

data_scaled <- scale(data, center = data_mean, scale = data_sd)

data

#-----------------------------------------------------------------------------------------------------------
#Lasso models
#-----------------------------------------------------------------------------------------------------------

y <- data[, 1]       # Extract the first column GDP
X <- data[, -1]      # Extract all columns except the first

#Training and test samples 
#70% of data used for training, 30% for testing. 
n <- nrow(X)
train_index <- floor(0.7 * n)
train <- 1:train_index
test <- (train_index + 1):n

#Lasso with timeslices
myTimeControl <- trainControl(method = "timeslice",
                              initialWindow = 16,  # 16 quarters (4 years)
                              horizon = 8,         # 8 quarters (2 years)
                              fixedWindow = TRUE ,
                              allowParallel = TRUE)


lasso_model <- train(X[train, ], y[train],
                     method = "glmnet",    # Lasso regression uses glmnet
                     trControl = myTimeControl,
                     metric='RMSE',
                     tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.05, by = 0.001))) 
lasso_model_pred <-  predict(lasso_model, X[test, ])

# Information Criteria
# BIC
lasso_bic_model <- ic.glmnet(X[train, ], y[train], crit = "bic")
lasso_pred_bic <- predict(lasso_bic_model, newdata = X[test, ])
# AIC
lasso_aic_model <- ic.glmnet(X[train, ], y[train], crit = "aic")
lasso_pred_aic <- predict(lasso_aic_model, newdata = X[test, ])


#See which variables were selected
# Variables selected for the Lasso time slice model
lasso_coeffs <- coef(lasso_model$finalModel, lasso_model$bestTune$lambda)

# Variables selected for BIC
coef.lasso.bic <- coef(lasso_bic_model)  # Get coefficients from the BIC-based model
selected.vars.bic <- which(coef.lasso.bic != 0)  # Non-zero coefficients
selected.vars.bic <- selected.vars.bic[names(selected.vars.bic) != "(Intercept)"]

# Variables selected for AIC
coef.lasso.aic <- coef(lasso_aic_model)  # Get coefficients from the AIC-based model
selected.vars.aic <- which(coef.lasso.aic != 0)  # Non-zero coefficients
selected.vars.aic <- selected.vars.aic[names(selected.vars.aic) != "(Intercept)"]

#-----------------------------------------------------------------------------------------------------------
#Evaluate predictions
#-----------------------------------------------------------------------------------------------------------

# Helper function to calculate RMSE and MAE
calculate_metrics <- function(actual, predicted) {
  rmse <- sqrt(mean((actual - predicted)^2))
  mae <- mean(abs(actual - predicted))
  list(RMSE = rmse, MAE = mae)
}

# Calculate metrics for each model
metrics_lasso_model <- calculate_metrics(y[test], lasso_model_pred)
metrics_lasso_bic <- calculate_metrics(y[test], lasso_pred_bic)
metrics_lasso_aic <- calculate_metrics(y[test], lasso_pred_aic)


# Count the number of variables selected for each model
num_vars_lasso_bic <- length(selected.vars.bic)
num_vars_lasso_aic <- length(selected.vars.aic)
num_vars_lasso_model <- sum(lasso_coeffs != 0) - 1  # Subtract 1 to exclude the intercept


# Calculate errors for each model
lasso_model_errors <- y[test] - lasso_model_pred
lasso_errors_bic <- y[test] - lasso_pred_bic
lasso_errors_aic <- y[test] - lasso_pred_aic

# Calculate q_t for each model
calculate_q_t <- function(errors, time_series) {
  n <- length(time_series)
  denominator <- (1 / (n - 1)) * sum(abs(diff(time_series)))
  q_t <- errors / denominator
  return(q_t)
}

# Calculate scaled errors (q_t) for each model
lasso_model_q_t <- calculate_q_t(lasso_model_errors, y[test])
lasso_q_t_bic <- calculate_q_t(lasso_errors_bic, y[test])
lasso_q_t_aic <- calculate_q_t(lasso_errors_aic, y[test])

# Calculate MASE for each model
MASE_lasso_model <- mean(abs(lasso_model_q_t))
MASE_lasso_bic <- mean(abs(lasso_q_t_bic))
MASE_lasso_aic <- mean(abs(lasso_q_t_aic))

# Retrieve the lambda for each model
lambda_lasso <- lasso_model$bestTune$lambda
lambda_bic <- lasso_bic_model[["lambda"]]
lambda_aic <- lasso_aic_model[["lambda"]]
  
# Combine the metrics into a data frame
results <- data.frame(
  Model = c("Lasso Timeslice", "Lasso BIC", "Lasso AIC"),
  RMSE = c(metrics_lasso_model$RMSE, metrics_lasso_bic$RMSE, metrics_lasso_aic$RMSE),
  MAE = c(metrics_lasso_model$MAE, metrics_lasso_bic$MAE, metrics_lasso_aic$MAE),
  MASE = c(MASE_lasso_model, MASE_lasso_bic, MASE_lasso_aic),
  Num_Variables = c(num_vars_lasso_model, num_vars_lasso_bic, num_vars_lasso_aic),
  Lambda = c(lambda_lasso, lambda_bic, lambda_aic)
)

print(results)

#Plot the actual vs forecasts of the different models

plot(time[test], y[test], type = "l", col = "blue", lwd = 2,  
     xlab = "Time", ylab = "Value",  
     xlim = range(time[test]), ylim = range(c(y[test], 
                                              lasso_model_pred, 
                                              lasso_pred_bic, 
                                              lasso_pred_aic)))

lines(time[test], lasso_model_pred, col = "orange", lwd = 1)
lines(time[test], lasso_pred_bic, col = "red", lwd = 2)
lines(time[test], lasso_pred_aic, col = "green", lwd = 1)
lines(time[test], y[test], col = rgb(0, 0, 1, alpha = 0.3), lwd = 1)


# Add a legend 
legend("bottomright", legend = c("Actual", "Lasso Timeslice", "Lasso BIC", "Lasso AIC"),
       col = c("blue", "orange", "red", "green"), 
       lty = 1, lwd = 2, bty = "n")  # cex controls the size of the legend text

