# sarimax_forecast.R
# Grid search for optimal SARIMAX parameters, fits model on training data,
# validates on held-out period, then forecasts streamflow under RCP 4.5 and 8.5.
#
# Best params from grid search: SARIMAX(6,1,1)(0,1,0)[365] — RMSE: 430.74 cfs
#
# Inputs: data/processed/combined_hist_vars.csv, combined_rcp45_vars.csv, combined_rcp85_vars.csv
# Outputs: data/processed/rcp45_streamflow_predictions.csv, rcp85_streamflow_predictions.csv

library(dplyr)
library(tseries)
library(lmtest)
library(stats)
library(forecast)
library(vars)
library(lubridate)

historic_vars <- read.csv("data/processed/combined_hist_vars.csv")
rcp45_vars <- read.csv("data/processed/combined_rcp45_vars.csv")
rcp85_vars <- read.csv("data/processed/combined_rcp85_vars.csv")

historic_vars$Date <- as.Date(historic_vars$Date)
rcp45_vars$Date <- as.Date(rcp45_vars$Date)
rcp85_vars$Date <- as.Date(rcp85_vars$Date)

# split into training and validation datasets
train_data <- historic_vars[historic_vars$Date <= as.Date("2018-12-31"), ]
val_data <- historic_vars[historic_vars$Date >= as.Date("2019-01-01"), ]

# define SARIMAX parameters for grid search
p_values <- 1:6
q_values <- 0:3
P_values <- 0:2
Q_values <- 0:2
d <- 1 # first order differencing for stationarity
D <- 1 # seasonal first-order differencing
S <- 365 # seasonal period

grid_search_results <- data.frame(
  p = integer(),
  d = integer(),
  q = integer(),
  P = integer(),
  D = integer(),
  Q = integer(),
  S = integer(),
  RMSE = numeric()
)

# grid search
for (p in p_values) {
  for (q in q_values) {
    for (P in P_values) {
      for (Q in Q_values) {
        cat("Testing params: p =",p,", q =",q,", P =",P,", Q =",Q,"\n")
        
        try({
          sarimax_model <- Arima(
            train_data$Streamflow,
            order = c(p, d, q),
            seasonal = c(P, D, Q, S),
            xreg = as.matrix(train_data[, c("Precip","SolRad","MinTemp","MaxTemp","VPD")])
          )
          
          val_predictions <- forecast(
            sarimax_model,
            h = nrow(val_data),
            xreg = as.matrix(val_data[, c("Precip","SolRad","MinTemp","MaxTemp","VPD")])
          )
          
          rmse_val <- sqrt(mean((val_data$Streamflow - val_predictions$mean)^2))
          
          grid_search_results <- rbind(
            grid_search_results,
            data.frame(p=p,d=d,q=q,P=P,D=D,Q=Q,S=S,RMSE=rmse_val)
          )
        }, silent = TRUE) # suppress errors from failed parameter combos
      }
    }
  }
}

best_params <- grid_search_results[which.min(grid_search_results$RMSE), ]
cat("Best Params:\n")
print(best_params)
#     p d q P D Q   S     RMSE
# 190 6 1 1 0 1 0 365 430.7403
grid_search_results <- grid_search_results[order(grid_search_results$RMSE), ]
write.csv(grid_search_results, "data/processed/sarimax_params_gridsearch_results.csv", row.names = FALSE)

# fit final model with best params on full historic data
p <- 6
d <- 1
q <- 1
P <- 0
D <- 1
Q <- 0
S <- 365

sarimax_model <- Arima(
  historic_vars$Streamflow,
  order = c(p, d, q),
  seasonal = c(P, D, Q, S),
  xreg = as.matrix(historic_vars[, c("Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")])
)

rcp45_combined <- rbind(
  historic_vars[, c("Date", "Streamflow", "Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")],
  cbind(Date = rcp45_vars$Date, Streamflow = NA, rcp45_vars[, c("Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")])
)

rcp45_forecast <- forecast(
  sarimax_model,
  h = nrow(rcp45_vars),
  xreg = as.matrix(rcp45_combined[nrow(historic_vars) + 1:nrow(rcp45_vars), c("Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")])
)

rcp45_results <- data.frame(
  Date = rcp45_combined$Date,
  Streamflow_Predicted = c(historic_vars$Streamflow, rcp45_forecast$mean)
)

write.csv(rcp45_results, "data/processed/rcp45_streamflow_predictions.csv", row.names = FALSE)

rcp85_combined <- rbind(
  historic_vars[, c("Date", "Streamflow", "Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")],
  cbind(Date = rcp85_vars$Date, Streamflow = NA, rcp85_vars[, c("Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")])
)

rcp85_forecast <- forecast(
  sarimax_model,
  h = nrow(rcp85_vars),
  xreg = as.matrix(rcp85_combined[nrow(historic_vars) + 1:nrow(rcp85_vars), c("Precip", "SolRad", "MinTemp", "MaxTemp", "VPD")])
)

rcp85_results <- data.frame(
  Date = rcp85_combined$Date,
  Streamflow_Predicted = c(historic_vars$Streamflow, rcp85_forecast$mean)
)

write.csv(rcp85_results, "data/processed/rcp85_streamflow_predictions.csv", row.names = FALSE)
