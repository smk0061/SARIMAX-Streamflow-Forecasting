# results_visualization.R
# Reads in SARIMAX predictions and generates comparison plots for
# historic actual vs modeled and future RCP projections.
#
# Inputs: data/processed/*_streamflow_predictions.csv

library(dplyr)
library(lubridate)
library(ggplot2)

historical_proj <- read.csv("data/processed/historic_streamflow_predictions.csv")
rcp45_proj <- read.csv("data/processed/rcp45_streamflow_predictions.csv")
rcp85_proj <- read.csv("data/processed/rcp85_streamflow_predictions.csv")

# aggregate to monthly for cleaner visualization
historical_proj <- historical_proj %>%
  mutate(YearMonth = floor_date(as.Date(Date), "month")) %>%
  group_by(YearMonth) %>%
  summarize(Streamflow_Actual = mean(Streamflow_Actual),
            Streamflow_Predicted = mean(Streamflow_Predicted))

ggplot(historical_proj, aes(x = YearMonth)) +
  geom_line(aes(y = Streamflow_Actual, color = "Actual")) +
  geom_line(aes(y = Streamflow_Predicted, color = "Predicted")) +
  labs(title = "Historic Streamflow vs Estimated Streamflow of Shavers Fork",
       x = "Year",
       y = "Streamflow (cfs)",
       color = "Legend") +
  scale_x_date(date_labels = "%Y", date_breaks = "5 years") +
  theme_minimal()

rcp45_proj <- rcp45_proj %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
         YearMonth = floor_date(Date, "month")) %>%
  filter(YearMonth >= as.Date("2024-01-01")) %>%
  group_by(YearMonth) %>%
  summarize(Streamflow_Predicted = mean(Streamflow_Predicted)) %>%
  mutate(Scenario = "RCP 4.5")

rcp85_proj <- rcp85_proj %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
         YearMonth = floor_date(Date, "month")) %>%
  filter(YearMonth >= as.Date("2024-01-01")) %>%
  group_by(YearMonth) %>%
  summarize(Streamflow_Predicted = mean(Streamflow_Predicted)) %>%
  mutate(Scenario = "RCP 8.5")

combined_proj <- bind_rows(rcp45_proj, rcp85_proj)

ggplot(combined_proj, aes(x = YearMonth, y = Streamflow_Predicted, color = Scenario)) +
  geom_line() +
  labs(title = "Projected Monthly Streamflow under RCP 4.5 and RCP 8.5 for Shavers Fork",
       x = "Year",
       y = "Streamflow (cfs)",
       color = "Scenario") +
  scale_x_date(date_labels = "%Y", date_breaks = "5 years") +
  theme_minimal()
