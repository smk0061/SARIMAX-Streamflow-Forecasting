# exploratory_analysis.R
# Exploratory data analysis: monthly averages, loess smoothing, leap day removal,
# stationarity testing (ADF/KPSS), STL decomposition, deseasonalization,
# standardized cross-variable comparison, PCA, and ACF/PACF analysis.
#
# This analysis informed the final model decisions:
#   - STL confirmed annual seasonality → S=365, D=1 in SARIMAX
#   - PCA showed variable redundancy → dropped humidity, wind, ET
#   - ACF/PACF identified lag structure → p=6 AR order
#   - Stationarity tests → differencing for precip and min temp
#
# Inputs: data/raw/ CSVs from data_extraction.R

library(dplyr)
library(tidyverse)
library(tidyplots)
library(ggplot2)
library(forecast)
library(lubridate)
library(tseries)
library(stats)
library(scales)
library(factoextra)
library(plotly)
library(lmtest)

# reading in historic daily data for Shavers Q and gridMET climate vars
shavers_Q <- read.csv("data/raw/shavers_daily_Q.csv")
precip_hist <- read.csv("data/raw/shavers_precipitation_amount.csv")
humidity_hist <- read.csv("data/raw/shavers_daily_mean_specific_humidity.csv")
solrad_hist <- read.csv("data/raw/shavers_daily_mean_shortwave_radiation_at_surface.csv")
mintemp_hist <- read.csv("data/raw/shavers_daily_minimum_temperature.csv")
maxtemp_hist <- read.csv("data/raw/shavers_daily_maximum_temperature.csv")
vpd_hist <- read.csv("data/raw/shavers_daily_mean_vapor_pressure_deficit.csv")

# reading in future daily data for MACA climate vars
precip_rcp45 <- read.csv("data/raw/shavers_pr_rcp45.csv")
precip_rcp85 <- read.csv("data/raw/shavers_pr_rcp85.csv")
humidity_rcp45 <- read.csv("data/raw/shavers_huss_rcp45.csv")
humidity_rcp85 <- read.csv("data/raw/shavers_huss_rcp85.csv")
solrad_rcp45 <- read.csv("data/raw/shavers_rsds_rcp45.csv")
solrad_rcp85 <- read.csv("data/raw/shavers_rsds_rcp85.csv")
mintemp_rcp45 <- read.csv("data/raw/shavers_tasmin_rcp45.csv")
mintemp_rcp85 <- read.csv("data/raw/shavers_tasmin_rcp85.csv")
maxtemp_rcp45 <- read.csv("data/raw/shavers_tasmax_rcp45.csv")
maxtemp_rcp85 <- read.csv("data/raw/shavers_tasmax_rcp85.csv")
vpd_rcp45 <- read.csv("data/raw/shavers_vpd_rcp45.csv")
vpd_rcp85 <- read.csv("data/raw/shavers_vpd_rcp85.csv")


# data prep
add_date_column <- function(data) {
  data %>% mutate(Date=ymd(paste(Year,Month,Day,sep="-")))
}

shavers_Q <- add_date_column(shavers_Q)
precip_hist <- add_date_column(precip_hist)
humidity_hist <- add_date_column(humidity_hist)
solrad_hist <- add_date_column(solrad_hist)
mintemp_hist <- add_date_column(mintemp_hist)
maxtemp_hist <- add_date_column(maxtemp_hist)
vpd_hist <- add_date_column(vpd_hist)
precip_rcp45 <- add_date_column(precip_rcp45)
precip_rcp85 <- add_date_column(precip_rcp85)
humidity_rcp45 <- add_date_column(humidity_rcp45)
humidity_rcp85 <- add_date_column(humidity_rcp85)
solrad_rcp45 <- add_date_column(solrad_rcp45)
solrad_rcp85 <- add_date_column(solrad_rcp85)
mintemp_rcp45 <- add_date_column(mintemp_rcp45)
mintemp_rcp85 <- add_date_column(mintemp_rcp85)
maxtemp_rcp45 <- add_date_column(maxtemp_rcp45)
maxtemp_rcp85 <- add_date_column(maxtemp_rcp85)
vpd_rcp45 <- add_date_column(vpd_rcp45)
vpd_rcp85 <- add_date_column(vpd_rcp85)


# monthly averages
calculate_monthly_stats <- function(data, value_col, stat_type = "mean") {
  data %>%
    group_by(Year, Month) %>%
    summarize(
      MonthlyStat = if (stat_type == "sum") sum(!!sym(value_col), na.rm = TRUE) else mean(!!sym(value_col), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Month) %>%
    summarize(
      MeanMonthlyStat = mean(MonthlyStat, na.rm = TRUE),
      .groups = "drop"
    )
}

q_monthly_avg <- calculate_monthly_stats(shavers_Q, "Q", "mean")
precip_monthly_totals <- calculate_monthly_stats(precip_hist, "Max", "sum")
humidity_monthly_avg <- calculate_monthly_stats(humidity_hist, "Min", "mean")
solar_monthly_avg <- calculate_monthly_stats(solrad_hist, "Max", "mean")
min_temp_monthly_avg <- calculate_monthly_stats(mintemp_hist, "Max", "mean")
max_temp_monthly_avg <- calculate_monthly_stats(maxtemp_hist, "Max", "mean")
vpd_monthly_avg <- calculate_monthly_stats(vpd_hist, "Min", "mean")

monthly_averages <- q_monthly_avg %>%
  rename(Q_cfs = MeanMonthlyStat) %>%
  left_join(precip_monthly_totals %>% rename(Precip_mm = MeanMonthlyStat), by = "Month") %>%
  left_join(humidity_monthly_avg %>% rename(Humidity_percent = MeanMonthlyStat), by = "Month") %>%
  left_join(solar_monthly_avg %>% rename(Solar_Wm2 = MeanMonthlyStat), by = "Month") %>%
  left_join(min_temp_monthly_avg %>% rename(Min_Temp_K = MeanMonthlyStat), by = "Month") %>%
  left_join(max_temp_monthly_avg %>% rename(Max_Temp_K = MeanMonthlyStat), by = "Month") %>%
  left_join(vpd_monthly_avg %>% rename(VPD_kPa = MeanMonthlyStat), by = "Month")

write.csv(monthly_averages, "data/raw/hist_monthly_averages_table.csv", row.names = FALSE)
print(monthly_averages)

# monthly average Q time series
q_monthly_avg <- shavers_Q %>%
  group_by(Year, Month) %>%
  summarize(MonthlyAvgQ = mean(Q, na.rm = TRUE), .groups = "drop") %>%
  mutate(Date = as.Date(paste(Year, Month, "01", sep = "-")))

ggplot(q_monthly_avg, aes(x = Date, y = MonthlyAvgQ)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Monthly Average Streamflow (Q) from 1/1/1998 to 12/31/2023",
       x = "Year", y = "Streamflow (cfs)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal()

# loess-smoothed Q
shavers_Q %>% 
  tidyplot(x=Date, y=Q) %>%
  add_curve_fit(method="loess", linewidth=1, se=FALSE) %>%
  add_title("Historical Discharge on Shavers Fork at Bowden") %>%
  adjust_x_axis(title="Date") %>%
  adjust_y_axis(title="Q (cfs)") %>%
  adjust_size(width=150,height=100)

# loess-smoothed climate variables
plot_climate_var <- function(data, y_label) {
  title <- paste("Mean", y_label, "Over Time\n(loess smoothing applied)")
  data %>%
    tidyplot(x = Date, y = Mean) %>%
    add_curve_fit(method = "loess", linewidth = 1, se=FALSE) %>%
    add_title(title) %>%
    adjust_x_axis(title = "Date") %>%
    adjust_y_axis(title = y_label) %>%
    adjust_size(width = 150, height = 100)
}

plot_climate_var(precip_hist, "Historical Precip")
plot_climate_var(precip_rcp45, "RCP 4.5 Precip")
plot_climate_var(precip_rcp85, "RCP 8.5 Precip")
plot_climate_var(humidity_hist, "Historical Humidity")
plot_climate_var(humidity_rcp45, "RCP 4.5 Humidity")
plot_climate_var(humidity_rcp85, "RCP 8.5 Humidity")
plot_climate_var(solrad_hist, "Historical Solar Radiation")
plot_climate_var(solrad_rcp45, "RCP 4.5 Solar Radiation")
plot_climate_var(solrad_rcp85, "RCP 8.5 Solar Radiation")
plot_climate_var(mintemp_hist, "Historical Min Temp")
plot_climate_var(mintemp_rcp45, "RCP 4.5 Min Temp")
plot_climate_var(mintemp_rcp85, "RCP 8.5 Min Temp")
plot_climate_var(maxtemp_hist, "Historical Max Temp")
plot_climate_var(maxtemp_rcp45, "RCP 4.5 Max Temp")
plot_climate_var(maxtemp_rcp85, "RCP 8.5 Max Temp")
plot_climate_var(vpd_hist, "Historical VPD")
plot_climate_var(vpd_rcp45, "RCP 4.5 VPD")
plot_climate_var(vpd_rcp85, "RCP 8.5 VPD")


## Remove Leap Days
# removing feb 29 so ts frequency can be set to 365

remove_leap_days <- function(data){
  data %>% filter(!(Month == 2 & Day == 29))
}

shavers_Q <- remove_leap_days(shavers_Q)
precip_hist <- remove_leap_days(precip_hist)
humidity_hist <- remove_leap_days(humidity_hist)
solrad_hist <- remove_leap_days(solrad_hist)
mintemp_hist <- remove_leap_days(mintemp_hist)
maxtemp_hist <- remove_leap_days(maxtemp_hist)
vpd_hist <- remove_leap_days(vpd_hist)
precip_rcp45 <- remove_leap_days(precip_rcp45)
precip_rcp85 <- remove_leap_days(precip_rcp85)
humidity_rcp45 <- remove_leap_days(humidity_rcp45)
humidity_rcp85 <- remove_leap_days(humidity_rcp85)
solrad_rcp45 <- remove_leap_days(solrad_rcp45)
solrad_rcp85 <- remove_leap_days(solrad_rcp85)
mintemp_rcp45 <- remove_leap_days(mintemp_rcp45)
mintemp_rcp85 <- remove_leap_days(mintemp_rcp85)
maxtemp_rcp45 <- remove_leap_days(maxtemp_rcp45)
maxtemp_rcp85 <- remove_leap_days(maxtemp_rcp85)
vpd_rcp45 <- remove_leap_days(vpd_rcp45)
vpd_rcp85 <- remove_leap_days(vpd_rcp85)


# stationarity Testing (ADF and KPSS)
shavers_Q.ts <- ts(shavers_Q$Q, start=c(1998,1), frequency=365)
precip_hist.ts <- ts(precip_hist$Mean, start=c(1998,1), frequency=365)
humidity_hist.ts <- ts(humidity_hist$Mean, start=c(1998,1), frequency=365)
solrad_hist.ts <- ts(solrad_hist$Mean, start=c(1998,1), frequency=365)
mintemp_hist.ts <- ts(mintemp_hist$Mean, start=c(1998,1), frequency=365)
maxtemp_hist.ts <- ts(maxtemp_hist$Mean, start=c(1998,1), frequency=365)
vpd_hist.ts <- ts(vpd_hist$Mean, start=c(1998,1), frequency=365)

precip_rcp45.ts <- ts(precip_rcp45$Mean, start=c(2024,1), frequency=365)
precip_rcp85.ts <- ts(precip_rcp85$Mean, start=c(2024,1), frequency=365)
humidity_rcp45.ts <- ts(humidity_rcp45$Mean, start=c(2024,1), frequency=365)
humidity_rcp85.ts <- ts(humidity_rcp85$Mean, start=c(2024,1), frequency=365)
solrad_rcp45.ts <- ts(solrad_rcp45$Mean, start=c(2024,1), frequency=365)
solrad_rcp85.ts <- ts(solrad_rcp85$Mean, start=c(2024,1), frequency=365)
mintemp_rcp45.ts <- ts(mintemp_rcp45$Mean, start=c(2024,1), frequency=365)
mintemp_rcp85.ts <- ts(mintemp_rcp85$Mean, start=c(2024,1), frequency=365)
maxtemp_rcp45.ts <- ts(maxtemp_rcp45$Mean, start=c(2024,1), frequency=365)
maxtemp_rcp85.ts <- ts(maxtemp_rcp85$Mean, start=c(2024,1), frequency=365)
vpd_rcp45.ts <- ts(vpd_rcp45$Mean, start=c(2024,1), frequency=365)
vpd_rcp85.ts <- ts(vpd_rcp85$Mean, start=c(2024,1), frequency=365)

stationarity_tests <- function(ts_data, ts_name) {
  adf_result <- adf.test(ts_data)
  adf_p_value <- adf_result$p.value
  adf_stationary <- ifelse(adf_p_value < 0.05, "Stationary", "Non-Stationary")
  
  kpss_result <- kpss.test(ts_data, null = "Trend")
  kpss_p_value <- kpss_result$p.value
  kpss_stationary <- ifelse(kpss_p_value < 0.05, "Non-Stationary", "Stationary")
  
  result_row <- data.frame(
    Variable = ts_name,
    ADF_Statistic = adf_result$statistic,
    ADF_p_value = adf_p_value,
    ADF_Stationary = adf_stationary,
    KPSS_Statistic = kpss_result$statistic,
    KPSS_p_value = kpss_p_value,
    KPSS_Stationary = kpss_stationary
  )
  return(result_row)
}

time_series_list <- list(
  shavers_Q.ts = shavers_Q.ts,
  precip_hist.ts = precip_hist.ts,
  humidity_hist.ts = humidity_hist.ts,
  solrad_hist.ts = solrad_hist.ts,
  mintemp_hist.ts = mintemp_hist.ts,
  maxtemp_hist.ts = maxtemp_hist.ts,
  vpd_hist.ts = vpd_hist.ts,
  precip_rcp45.ts = precip_rcp45.ts,
  precip_rcp85.ts = precip_rcp85.ts,
  humidity_rcp45.ts = humidity_rcp45.ts,
  humidity_rcp85.ts = humidity_rcp85.ts,
  solrad_rcp45.ts = solrad_rcp45.ts,
  solrad_rcp85.ts = solrad_rcp85.ts,
  mintemp_rcp45.ts = mintemp_rcp45.ts,
  mintemp_rcp85.ts = mintemp_rcp85.ts,
  maxtemp_rcp45.ts = maxtemp_rcp45.ts,
  maxtemp_rcp85.ts = maxtemp_rcp85.ts,
  vpd_rcp45.ts = vpd_rcp45.ts,
  vpd_rcp85.ts = vpd_rcp85.ts
)

stationarity_results <- data.frame()
for (ts_name in names(time_series_list)) {
  ts_data <- time_series_list[[ts_name]]
  result_row <- stationarity_tests(ts_data, ts_name)
  stationarity_results <- rbind(stationarity_results, result_row)
}
print(stationarity_results)


## STL Decomposition

stl_decompose_and_plot <- function(ts_data, main_title_prefix) {
  stl_decomp <- stl(ts_data, s.window = "periodic")
  plot(stl_decomp, main = paste("STL Decomposition of", main_title_prefix))
  return(stl_decomp)
}

stl_decompose_and_plot(shavers_Q.ts, "Discharge (Historical)")
stl_decompose_and_plot(precip_hist.ts, "Precipitation (Historical)")
stl_decompose_and_plot(precip_rcp45.ts, "Precipitation (RCP 4.5)")
stl_decompose_and_plot(precip_rcp85.ts, "Precipitation (RCP 8.5)")
stl_decompose_and_plot(humidity_hist.ts, "Specific Humidity (Historical)")
stl_decompose_and_plot(humidity_rcp45.ts, "Specific Humidity (RCP 4.5)")
stl_decompose_and_plot(humidity_rcp85.ts, "Specific Humidity (RCP 8.5)")
stl_decompose_and_plot(solrad_hist.ts, "Solar Radiation (Historical)")
stl_decompose_and_plot(solrad_rcp45.ts, "Solar Radiation (RCP 4.5)")
stl_decompose_and_plot(solrad_rcp85.ts, "Solar Radiation (RCP 8.5)")
stl_decompose_and_plot(mintemp_hist.ts, "Min Temperature (Historical)")
stl_decompose_and_plot(mintemp_rcp45.ts, "Min Temperature (RCP 4.5)")
stl_decompose_and_plot(mintemp_rcp85.ts, "Min Temperature (RCP 8.5)")
stl_decompose_and_plot(maxtemp_hist.ts, "Max Temperature (Historical)")
stl_decompose_and_plot(maxtemp_rcp45.ts, "Max Temperature (RCP 4.5)")
stl_decompose_and_plot(maxtemp_rcp85.ts, "Max Temperature (RCP 8.5)")
stl_decompose_and_plot(vpd_hist.ts, "VPD (Historical)")
stl_decompose_and_plot(vpd_rcp45.ts, "VPD (RCP 4.5)")
stl_decompose_and_plot(vpd_rcp85.ts, "VPD (RCP 8.5)")


## Deseasonalization
# remove seasonal component so only standardizing trend and residuals

deseasonalize_var <- function(ts_data) {
  stl_decomp <- stl(ts_data, s.window = "periodic")
  deseasonalize_data <- stl_decomp$time.series[, "trend"] + stl_decomp$time.series[, "remainder"]
  return(deseasonalize_data)
}

shavers_Q$deseasonalize_mean <- deseasonalize_var(shavers_Q.ts)
precip_hist$deseasonalize_mean <- deseasonalize_var(precip_hist.ts)
precip_rcp45$deseasonalize_mean <- deseasonalize_var(precip_rcp45.ts)
precip_rcp85$deseasonalize_mean <- deseasonalize_var(precip_rcp85.ts)
humidity_hist$deseasonalized_mean <- deseasonalize_var(humidity_hist.ts)
humidity_rcp45$deseasonalized_mean <- deseasonalize_var(humidity_rcp45.ts)
humidity_rcp85$deseasonalized_mean <- deseasonalize_var(humidity_rcp85.ts)
solrad_hist$deseasonalized_mean <- deseasonalize_var(solrad_hist.ts)
solrad_rcp45$deseasonalized_mean <- deseasonalize_var(solrad_rcp45.ts)
solrad_rcp85$deseasonalized_mean <- deseasonalize_var(solrad_rcp85.ts)
mintemp_hist$deseasonalized_mean <- deseasonalize_var(mintemp_hist.ts)
mintemp_rcp45$deseasonalized_mean <- deseasonalize_var(mintemp_rcp45.ts)
mintemp_rcp85$deseasonalized_mean <- deseasonalize_var(mintemp_rcp85.ts)
maxtemp_hist$deseasonalized_mean <- deseasonalize_var(maxtemp_hist.ts)
maxtemp_rcp45$deseasonalized_mean <- deseasonalize_var(maxtemp_rcp45.ts)
maxtemp_rcp85$deseasonalized_mean <- deseasonalize_var(maxtemp_rcp85.ts)
vpd_hist$deseasonalized_mean <- deseasonalize_var(vpd_hist.ts)
vpd_rcp45$deseasonalized_mean <- deseasonalize_var(vpd_rcp45.ts)
vpd_rcp85$deseasonalized_mean <- deseasonalize_var(vpd_rcp85.ts)

# deseasonalized plots
shavers_Q %>%
  tidyplot(x = Date, y = deseasonalize_mean, geom = "line") %>%
  add_line() %>%
  add_title("Deseasonalized Discharge (Historical)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)

precip_hist %>%
  tidyplot(x = Date, y = deseasonalize_mean, geom = "line") %>%
  add_line(linewidth=0.01) %>%
  add_title("Deseasonalized Precipitation (Historical)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)

precip_rcp45 %>%
  tidyplot(x = Date, y = deseasonalize_mean, geom = "line") %>%
  add_line(linewidth=0.01) %>%
  add_title("Deseasonalized Precipitation (RCP 4.5)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)

precip_rcp85 %>%
  tidyplot(x = Date, y = deseasonalize_mean, geom = "line") %>%
  add_line(linewidth=0.01) %>%
  add_title("Deseasonalized Precipitation (RCP 8.5)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)

humidity_hist %>%
  tidyplot(x = Date, y = deseasonalized_mean, geom = "line") %>%
  add_line(linewidth=0.01) %>%
  add_title("Deseasonalized Specific Humidity (Historical)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)

solrad_hist %>%
  tidyplot(x = Date, y = deseasonalized_mean, geom = "line") %>%
  add_line(linewidth=0.01) %>%
  add_title("Deseasonalized Solar Radiation (Historical)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)

mintemp_rcp85 %>%
  tidyplot(x = Date, y = deseasonalized_mean, geom = "line") %>%
  add_line(linewidth=0.01) %>%
  add_title("Deseasonalized Min Temperature (RCP 8.5)") %>%
  adjust_x_axis(title = "Date") %>%
  adjust_y_axis(title = "Deseasonalized Mean") %>%
  adjust_size(width = 225, height = 125)


# standardized cross var comparison
# read in all historic climate as a list for standardization
daily_Q <- read.csv("data/raw/shavers_daily_Q.csv")
daily_precip <- read.csv("data/raw/shavers_precipitation_amount.csv")
daily_humidity <- read.csv("data/raw/shavers_daily_mean_specific_humidity.csv")
daily_solrad <- read.csv("data/raw/shavers_daily_mean_shortwave_radiation_at_surface.csv")
daily_mintemp <- read.csv("data/raw/shavers_daily_minimum_temperature.csv")
daily_maxtemp <- read.csv("data/raw/shavers_daily_maximum_temperature.csv")
daily_windspd <- read.csv("data/raw/shavers_daily_mean_wind_speed.csv")
daily_evapotrans <- read.csv("data/raw/shavers_daily_mean_reference_evapotranspiration_grass.csv")
daily_vpd <- read.csv("data/raw/shavers_daily_mean_vapor_pressure_deficit.csv")

climate_data_list <- list(
  daily_precip = daily_precip,
  daily_humidity = daily_humidity,
  daily_solrad = daily_solrad,
  daily_mintemp = daily_mintemp,
  daily_maxtemp = daily_maxtemp,
  daily_windspd = daily_windspd,
  daily_evapotrans = daily_evapotrans,
  daily_vpd = daily_vpd
)

# standardize all climate vars and Q
standardize_ts <- function(data, value_col) {
  data %>%
    mutate(Standardized = (get(value_col) - mean(get(value_col), na.rm = TRUE)) / sd(get(value_col), na.rm = TRUE))
}

standard_mean_var <- data.frame()
for(var_name in names(climate_data_list)){
  climate_data <- climate_data_list[[var_name]]
  standard_data <- standardize_ts(climate_data, "Mean")
  standard_data$Variable <- var_name
  standard_data <- standard_data %>% select(Year, Month, Day, Variable, Standardized)
  standard_mean_var <- rbind(standard_mean_var, standard_data)
}

standard_Q <- standardize_ts(daily_Q, "Q")
standard_Q$Variable <- "Streamflow"
standard_Q <- standard_Q %>% select(Year, Month, Day, Variable, Standardized)

standard_combined <- rbind(standard_mean_var, standard_Q)
standard_combined$Date <- as.Date(with(standard_combined, paste(Year, Month, Day, sep = "-")), "%Y-%m-%d")

# smoothed standardized comparison plot
ggplot(standard_combined, aes(x = Date, y = Standardized, color = Variable)) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, size = 1) +
  scale_color_brewer(palette = "Paired") +
  labs(
    title = "Smoothed Standardized Climate Variables and Streamflow for Shavers Fork Watershed (1998-2023)",
    subtitle = "LOESS smoothing applied to standardized daily means and streamflow for easier trend comparison",
    x = "Date",
    y = "Standardized Value",
    color = "Variable"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 10)
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "5 years") +
  scale_y_continuous(labels = comma)

# monthly means of standardized data
standard_monthly_avg <- standard_combined %>%
  group_by(Year, Month, Variable) %>%
  summarise(Monthly_Mean_Standardized = mean(Standardized, na.rm = TRUE)) %>%
  ungroup()

standard_monthly_avg$Date <- as.Date(paste(standard_monthly_avg$Year, standard_monthly_avg$Month, "01", sep = "-"), format = "%Y-%m-%d")

ggplot(standard_monthly_avg, aes(x = Date, y = Monthly_Mean_Standardized, color = Variable)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Monthly Mean Standardized Climate Variables and Streamflow",
    subtitle = "Standardized Monthly Means for Shavers Fork Watershed (1998-2023)",
    x = "Year-Month",
    y = "Monthly Mean Standardized Value",
    color = "Variable"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  scale_y_continuous(labels = scales::comma)


# PCA
# pivot standardized data wide for PCA
standard_combined_pca_ready <- standard_combined %>%
  select(Year, Month, Day, Variable, Standardized) %>%
  pivot_wider(names_from = Variable, values_from = Standardized) %>%
  select(-Year, -Month, -Day)

pca_result <- prcomp(standard_combined_pca_ready, center = TRUE, scale. = TRUE)
summary(pca_result)

# scree plot
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# PCA loadings
pca_loadings <- pca_result$rotation
print(pca_loadings)

# biplot
fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "darkviolet",
                col.ind = "black",
                alpha.ind = 0.5,
                geom.ind = "point"
)

pca_scores <- as.data.frame(pca_result$x)

# 3d PCA plot
plot_ly(data = pca_scores, 
        x = ~PC1, y = ~PC2, z = ~PC3, 
        type = 'scatter3d', mode = 'markers',
        marker = list(size = 3, color = ~PC1, colorscale = "Viridis")) %>%
  layout(scene = list(
    xaxis = list(title = 'PC1'),
    yaxis = list(title = 'PC2'),
    zaxis = list(title = 'PC3')
  ),
  title = "3D PCA Plot of Climate Variables and Streamflow")

# parallel coordinates plot
pca_scores_long <- pca_scores %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "Component", values_to = "Value")

ggplot(pca_scores_long, aes(x = Component, y = Value, group = ID, color = Component)) +
  geom_line(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Parallel Coordinates Plot of First Nine Principal Components",
       x = "Principal Components", y = "Value") +
  theme(legend.position = "none") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "yellow")

print(pca_loadings)


## ACF and PACF Analysis

variables <- c("Streamflow", names(climate_data_list))

lag_periods <- list(
  Annual = 365,
  QuarterAnnual = 91,
  Monthly = 30,
  TenDays = 10,
  FiveDays = 5,
  ThreeDays = 3
)

for (var_name in variables) {
  if (var_name == "Streamflow") {
    data <- daily_Q
    value_column <- "Q"
  } else {
    data <- climate_data_list[[var_name]]
    value_column <- "Mean"
  }
  
  values <- data[[value_column]]
  
  for (lag_name in names(lag_periods)) {
    lag_max <- lag_periods[[lag_name]]
    
    cat("\nACF and PACF for", var_name, "using lag period:", lag_name, "\n")
    
    par(mfrow = c(1, 2))
    acf(values, lag.max = lag_max, main = paste("ACF -", var_name, "-", lag_name))
    pacf(values, lag.max = lag_max, main = paste("PACF -", var_name, "-", lag_name))
  }
}

# significant lag detection
for (var_name in variables) {
  if (var_name == "Streamflow") {
    data <- daily_Q
    value_column <- "Q"
  } else {
    data <- climate_data_list[[var_name]]
    value_column <- "Mean"
  }
  
  values <- data[[value_column]]
  N <- length(values)
  threshold <- 1.96 / sqrt(N)
  
  for (lag_name in names(lag_periods)) {
    lag_max <- lag_periods[[lag_name]]
    
    acf_result <- acf(values, lag.max = lag_max, plot = FALSE)
    pacf_result <- pacf(values, lag.max = lag_max, plot = FALSE)
    
    significant_acf_lags <- which(abs(acf_result$acf) > threshold)
    significant_pacf_lags <- which(abs(pacf_result$acf) > threshold)
    
    if (length(significant_acf_lags) > 0) {
      cat("significant ACF lags for", var_name, "(", lag_name, "):", significant_acf_lags - 1, "\n")
    }
    if (length(significant_pacf_lags) > 0) {
      cat("significant PACF lags for", var_name, "(", lag_name, "):", significant_pacf_lags - 1, "\n")
    }
  }
}
