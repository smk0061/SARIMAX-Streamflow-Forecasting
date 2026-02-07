# data_preparation.R
# Cleans extracted climate data, selects appropriate zonal statistics,
# runs stationarity tests, applies differencing where needed, runs Granger
# causality, standardizes using historic parameters (applied to RCP data),
# and combines into model-ready CSVs
#
# Inputs: data/raw/ CSVs from data_extraction.R
# Outputs: data/processed/combined_hist_vars.csv, combined_rcp45_vars.csv, combined_rcp85_vars.csv

library(dplyr)
library(tseries)
library(lmtest)
library(stats)
library(forecast)
library(vars)
library(lubridate)
library(ggplot2)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# historic climate
shavers_Q <- read.csv("data/raw/shavers_daily_Q.csv")
precip_hist <- read.csv("data/raw/shavers_precipitation_amount.csv")
solrad_hist <- read.csv("data/raw/shavers_daily_mean_shortwave_radiation_at_surface.csv")
mintemp_hist <- read.csv("data/raw/shavers_daily_minimum_temperature.csv")
maxtemp_hist <- read.csv("data/raw/shavers_daily_maximum_temperature.csv")
vpd_hist <- read.csv("data/raw/shavers_daily_mean_vapor_pressure_deficit.csv")

# future climate (MACA)
precip_rcp45 <- read.csv("data/raw/shavers_pr_rcp45.csv")
precip_rcp85 <- read.csv("data/raw/shavers_pr_rcp85.csv")
solrad_rcp45 <- read.csv("data/raw/shavers_rsds_rcp45.csv")
solrad_rcp85 <- read.csv("data/raw/shavers_rsds_rcp85.csv")
mintemp_rcp45 <- read.csv("data/raw/shavers_tasmin_rcp45.csv")
mintemp_rcp85 <- read.csv("data/raw/shavers_tasmin_rcp85.csv")
maxtemp_rcp45 <- read.csv("data/raw/shavers_tasmax_rcp45.csv")
maxtemp_rcp85 <- read.csv("data/raw/shavers_tasmax_rcp85.csv")
vpd_rcp45 <- read.csv("data/raw/shavers_vpd_rcp45.csv")
vpd_rcp85 <- read.csv("data/raw/shavers_vpd_rcp85.csv")

add_date_column <- function(data) {
  data %>% mutate(Date=ymd(paste(Year,Month,Day,sep="-")))
}

shavers_Q <- add_date_column(shavers_Q)
precip_hist <- add_date_column(precip_hist)
solrad_hist <- add_date_column(solrad_hist)
mintemp_hist <- add_date_column(mintemp_hist)
maxtemp_hist <- add_date_column(maxtemp_hist)
vpd_hist <- add_date_column(vpd_hist)
precip_rcp45 <- add_date_column(precip_rcp45)
precip_rcp85 <- add_date_column(precip_rcp85)
solrad_rcp45 <- add_date_column(solrad_rcp45)
solrad_rcp85 <- add_date_column(solrad_rcp85)
mintemp_rcp45 <- add_date_column(mintemp_rcp45)
mintemp_rcp85 <- add_date_column(mintemp_rcp85)
maxtemp_rcp45 <- add_date_column(maxtemp_rcp45)
maxtemp_rcp85 <- add_date_column(maxtemp_rcp85)
vpd_rcp45 <- add_date_column(vpd_rcp45)
vpd_rcp85 <- add_date_column(vpd_rcp85)

# select appropriate zonal statistic per variable
# Q = Q, Precip = Max, Sol rad = Max, Min temp = Max, Max temp = Max, VPD = Min

clean_and_save <- function(data, value_col, output_file) {
  cleaned_data <- data %>%
    select(Year, Month, Day, Date, all_of(value_col)) %>%
    rename(Value = all_of(value_col))
  
  write.csv(cleaned_data, output_file, row.names = FALSE)
}

clean_and_save(shavers_Q, "Q", "data/processed/cleaned_shavers_Q.csv")
clean_and_save(precip_hist, "Max", "data/processed/cleaned_precip_hist.csv")
clean_and_save(solrad_hist, "Max", "data/processed/cleaned_solrad_hist.csv")
clean_and_save(mintemp_hist, "Max", "data/processed/cleaned_mintemp_hist.csv")
clean_and_save(maxtemp_hist, "Max", "data/processed/cleaned_maxtemp_hist.csv")
clean_and_save(vpd_hist, "Min", "data/processed/cleaned_vpd_hist.csv")
clean_and_save(precip_rcp45, "Max", "data/processed/cleaned_precip_rcp45.csv")
clean_and_save(precip_rcp85, "Max", "data/processed/cleaned_precip_rcp85.csv")
clean_and_save(solrad_rcp45, "Max", "data/processed/cleaned_solrad_rcp45.csv")
clean_and_save(solrad_rcp85, "Max", "data/processed/cleaned_solrad_rcp85.csv")
clean_and_save(mintemp_rcp45, "Max", "data/processed/cleaned_mintemp_rcp45.csv")
clean_and_save(mintemp_rcp85, "Max", "data/processed/cleaned_mintemp_rcp85.csv")
clean_and_save(maxtemp_rcp45, "Max", "data/processed/cleaned_maxtemp_rcp45.csv")
clean_and_save(maxtemp_rcp85, "Max", "data/processed/cleaned_maxtemp_rcp85.csv")
clean_and_save(vpd_rcp45, "Min", "data/processed/cleaned_vpd_rcp45.csv")
clean_and_save(vpd_rcp85, "Min", "data/processed/cleaned_vpd_rcp85.csv")

# read in cleaned historical data
shavers_Q <- read.csv("data/processed/cleaned_shavers_Q.csv")
precip_hist <- read.csv("data/processed/cleaned_precip_hist.csv")
solrad_hist <- read.csv("data/processed/cleaned_solrad_hist.csv")
mintemp_hist <- read.csv("data/processed/cleaned_mintemp_hist.csv")
maxtemp_hist <- read.csv("data/processed/cleaned_maxtemp_hist.csv")
vpd_hist <- read.csv("data/processed/cleaned_vpd_hist.csv")

# ADF and KPSS tests on cleaned data
test_stationarity <- function(data, value_col, variable_name) {
  series <- data[[value_col]]
  cat("\nresults for:", variable_name, "\n")
  
  adf_result <- adf.test(series, alternative = "stationary")
  cat("ADF test:\n")
  print(adf_result)
  
  kpss_result <- kpss.test(series, null = "Level")
  cat("\nKPSS level test:\n")
  print(kpss_result)
  
  kpss_result2 <- kpss.test(series, null = "Trend")
  cat("\nKPSS trend test:\n")
  print(kpss_result2)
}

test_stationarity(shavers_Q, "Value", "Streamflow (Q)")
test_stationarity(precip_hist, "Value", "Precipitation")
test_stationarity(solrad_hist, "Value", "Solar Radiation")
test_stationarity(mintemp_hist, "Value", "Minimum Temperature")
test_stationarity(maxtemp_hist, "Value", "Maximum Temperature")
test_stationarity(vpd_hist, "Value", "Vapor Pressure Deficit")

# differencing precip and min temp for non-stationarity
apply_differencing <- function(data, value_col, output_file) {
  data <- data %>% mutate(DifferencedValue = c(NA, diff(!!sym(value_col))))
  data <- data %>% filter(!is.na(DifferencedValue))
  write.csv(data, output_file, row.names = FALSE)
}

apply_differencing(precip_hist, "Value", "data/processed/differenced_precipitation.csv")
apply_differencing(mintemp_hist, "Value", "data/processed/differenced_minimum_temperature.csv")

precip_hist_diff <- read.csv("data/processed/differenced_precipitation.csv")
mintemp_hist_diff<- read.csv("data/processed/differenced_minimum_temperature.csv")

test_stationarity(precip_hist_diff, "DifferencedValue", "Differenced Precipitation")
test_stationarity(mintemp_hist_diff, "DifferencedValue", "Differenced Minimum Temperature")

# Remove the first date (1/1/1998) from all datasets because of differencing
align_datasets <- function(data, date_col) {
  data %>% filter(!!sym(date_col) != "1998-01-01")
}

shavers_Q2 <- align_datasets(shavers_Q, "Date")
precip_hist2 <- align_datasets(precip_hist, "Date")
solrad_hist2 <- align_datasets(solrad_hist, "Date")
mintemp_hist2 <- align_datasets(mintemp_hist, "Date")
maxtemp_hist2 <- align_datasets(maxtemp_hist, "Date")
vpd_hist2 <- align_datasets(vpd_hist, "Date")
precip_hist_diff <- align_datasets(precip_hist_diff, "Date")
mintemp_hist_diff <- align_datasets(mintemp_hist_diff, "Date")

# combine vars for Granger test
combined_hist_data <- data.frame(
  Date = shavers_Q2$Date,
  Q = shavers_Q2$Value,
  Precip = precip_hist_diff$DifferencedValue,
  SolarRad = solrad_hist2$Value,
  MinTemp = mintemp_hist_diff$DifferencedValue,
  MaxTemp = maxtemp_hist2$Value,
  VPD = vpd_hist2$Value
)

# ACF and PACF for streamflow
acf_pacf_plot <- function(data, value_col, variable_name) {
  series <- data[[value_col]]
  
  par(mfrow = c(1, 2))
  acf(series, main = paste("ACF of", variable_name))
  pacf(series, main = paste("PACF of", variable_name))
  
  return(series)
}

streamflow_series <- acf_pacf_plot(shavers_Q2, "Value", "Streamflow (Q)")

# Granger causality testing across lags 1-30
test_granger_lags <- function(data, max_lag = 30) {
  climate_vars <- c("Precip", "SolarRad", "MinTemp", "MaxTemp", "VPD")
  
  for(var in climate_vars) {
    cat("\ngranger causality results for", var, "\n")
    cat("testing lags 1 to", max_lag, "\n")
    
    results <- data.frame(
      lag = 1:max_lag,
      p_value = NA
    )
    
    for(i in 1:max_lag) {
      g_test <- grangertest(Q ~ get(var), order = i, data = data)
      results$p_value[i] <- g_test$`Pr(>F)`[2]
    }
    
    sig_lags <- results[results$p_value < 0.05, ]
    cat("\nsignificant lags (p < 0.05):\n")
    print(sig_lags)
  }
}

test_granger_lags(combined_hist_data)

# read in future RCP data
precip_45 <- read.csv("data/processed/cleaned_precip_rcp45.csv")
precip_85 <- read.csv("data/processed/cleaned_precip_rcp85.csv")
solrad_45 <- read.csv("data/processed/cleaned_solrad_rcp45.csv")
solrad_85 <- read.csv("data/processed/cleaned_solrad_rcp85.csv")
mintemp_45 <- read.csv("data/processed/cleaned_mintemp_rcp45.csv")
mintemp_85 <- read.csv("data/processed/cleaned_mintemp_rcp85.csv")
maxtemp_45 <- read.csv("data/processed/cleaned_maxtemp_rcp45.csv")
maxtemp_85 <- read.csv("data/processed/cleaned_maxtemp_rcp85.csv")
vpd_45 <- read.csv("data/processed/cleaned_vpd_rcp45.csv")
vpd_85 <- read.csv("data/processed/cleaned_vpd_rcp85.csv")

# standardize historic data, save mean/sd for applying to RCP data
standardize <- function(data, value_col) {
  mean_val <- mean(data[[value_col]], na.rm = TRUE)
  sd_val <- sd(data[[value_col]], na.rm = TRUE)
  data <- data %>% mutate(Standardized = (.[[value_col]] - mean_val) / sd_val)
  return(list(data = data, mean = mean_val, sd = sd_val))
}

precip_hist_std <- standardize(precip_hist, "Value")
solrad_hist_std <- standardize(solrad_hist, "Value")
mintemp_hist_std <- standardize(mintemp_hist, "Value")
maxtemp_hist_std <- standardize(maxtemp_hist, "Value")
vpd_hist_std <- standardize(vpd_hist, "Value")

# apply historic standardization parameters to RCP data
apply_standardization <- function(data, mean_val, sd_val, value_col) {
  data <- data %>% mutate(Standardized = (.[[value_col]] - mean_val) / sd_val)
  return(data)
}

precip_rcp45_std <- apply_standardization(precip_45, precip_hist_std$mean, precip_hist_std$sd, "Value")
precip_rcp85_std <- apply_standardization(precip_85, precip_hist_std$mean, precip_hist_std$sd, "Value")
solrad_rcp45_std <- apply_standardization(solrad_45, solrad_hist_std$mean, solrad_hist_std$sd, "Value")
solrad_rcp85_std <- apply_standardization(solrad_85, solrad_hist_std$mean, solrad_hist_std$sd, "Value")
mintemp_rcp45_std <- apply_standardization(mintemp_45, mintemp_hist_std$mean, mintemp_hist_std$sd, "Value")
mintemp_rcp85_std <- apply_standardization(mintemp_85, mintemp_hist_std$mean, mintemp_hist_std$sd, "Value")
maxtemp_rcp45_std <- apply_standardization(maxtemp_45, maxtemp_hist_std$mean, maxtemp_hist_std$sd, "Value")
maxtemp_rcp85_std <- apply_standardization(maxtemp_85, maxtemp_hist_std$mean, maxtemp_hist_std$sd, "Value")
vpd_rcp45_std <- apply_standardization(vpd_45, vpd_hist_std$mean, vpd_hist_std$sd, "Value")
vpd_rcp85_std <- apply_standardization(vpd_85, vpd_hist_std$mean, vpd_hist_std$sd, "Value")

combined_hist_vars <- data.frame(
  Date = shavers_Q$Date,
  Streamflow = shavers_Q$Value,
  Precip = precip_hist_std$data$Standardized,
  SolRad = solrad_hist_std$data$Standardized,
  MinTemp = mintemp_hist_std$data$Standardized,
  MaxTemp = maxtemp_hist_std$data$Standardized,
  VPD = vpd_hist_std$data$Standardized
)

combined_rcp45_vars <- data.frame(
  Date = precip_45$Date,
  Precip = precip_rcp45_std$Standardized,
  SolRad = solrad_rcp45_std$Standardized,
  MinTemp = mintemp_rcp45_std$Standardized,
  MaxTemp = maxtemp_rcp45_std$Standardized,
  VPD = vpd_rcp45_std$Standardized
)

combined_rcp85_vars <- data.frame(
  Date = precip_85$Date,
  Precip = precip_rcp85_std$Standardized,
  SolRad = solrad_rcp85_std$Standardized,
  MinTemp = mintemp_rcp85_std$Standardized,
  MaxTemp = maxtemp_rcp85_std$Standardized,
  VPD = vpd_rcp85_std$Standardized
)

write.csv(combined_hist_vars, "data/processed/combined_hist_vars.csv", row.names = FALSE)
write.csv(combined_rcp45_vars, "data/processed/combined_rcp45_vars.csv", row.names = FALSE)
write.csv(combined_rcp85_vars, "data/processed/combined_rcp85_vars.csv", row.names = FALSE)
