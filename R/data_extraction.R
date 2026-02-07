# data_extraction.R
# Extracts USGS streamflow, gridMET historic climate, and MACA future climate 
# data for the Shavers Fork HUC-10 watershed. Clips rasters to watershed boundary
# and calculates daily zonal statistics. Exports to data/raw/
#
# Data sources:
#   USGS NWIS streamflow: https://waterdata.usgs.gov/nwis
#   gridMET climate: https://www.climatologylab.org/gridmet.html
#   MACA projections: https://www.climatologylab.org/maca.html
#
# Note: it will take several hours to donwload MACA data.

library(devtools)
library(dplyr)
library(lubridate)
library(dataRetrieval)
# devtools::install_github("mikejohnson51/climateR")
library(climateR)
library(terra)
library(raster)
library(sf)
library(sp)
library(tidyr)
library(ggplot2)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# shavers fork at bowden stream gauge
shavers_low <- "03068800"
startDate <- "1998-01-01"
endDate <- "2023-12-31"

# read in USGS streamflow data ("00060" = Q)
shavers_daily_Q <- readNWISdv(siteNumber = shavers_low, parameterCd = "00060", startDate = startDate, endDate = endDate)

# aggregate into monthly averages
shavers_monthly_Q <- shavers_daily_Q %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarise(Monthly_Discharge = mean(X_00060_00003, na.rm = TRUE)) %>%
  ungroup()

# organize daily Q data
shavers_daily_Q <- shavers_daily_Q %>%
  mutate(Year = year(Date), Month = month(Date), Day = day(Date), Q = X_00060_00003) %>%
  select(Year, Month, Day, Q)

write.csv(shavers_daily_Q, "data/raw/shavers_daily_Q.csv", row.names=FALSE)

ggplot(shavers_monthly_Q, aes(x = as.Date(paste(Year, Month, "01", sep = "-")), y = Monthly_Discharge)) +
  geom_line(color = "blue", group = 1)+
  labs(title = "Monthly Average Q for Shavers Fork at Bowden (03068800)",
       x = "Date", y = "Discharge (cfs)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_minimal()

ggplot(shavers_monthly_Q, aes(x = factor(Month), y = Monthly_Discharge))+
  geom_boxplot(fill = "lightblue") +
  labs(title = "Seasonal Trends in Monthly Discharge (Downstream Gauge)",
       x = "Month",y = "Discharge (cfs)") +
  theme_minimal()

# stream gauge coordinates
shavers_lower_gauge_lat <- 38.91316335
shavers_lower_gauge_long <- -79.770342

# HUC 10 Shavers Fork watershed shapefile
shavers_shp <- st_read("shp/shaversfork_huc10.shp")
shavers_prj <- st_transform(shavers_shp, crs = 4326)
plot(st_geometry(shavers_prj))

gauge_point <- SpatialPoints(cbind(shavers_lower_gauge_long, shavers_lower_gauge_lat))
points(gauge_point, col = "red", pch = 19, cex = 1.5)
text(shavers_lower_gauge_long, shavers_lower_gauge_lat, labels = "USGS: 03068800", pos = 4, col = "red")

# gridMET variables:
# pr (Precipitation), tmmx (Max Temp), tmmn (Min Temp), vs (Wind Speed),
# sph (Specific Humidity), srad (Solar Radiation), vpd (VPD), pet (Evapotranspiration)
gridMET_daily_vars <- getGridMET(
  AOI = shavers_prj,
  varname = c("pr", "tmmx", "tmmn", "vs", "sph", "srad", "vpd", "pet"),
  startDate = startDate,
  endDate = endDate
)

# write to 32-bit float temp files to manage memory
gridMET_daily_vars <- lapply(gridMET_daily_vars, function(raster_layer) {
  temp_file <- tempfile(fileext = ".tif")
  raster_layer <- terra::writeRaster(raster_layer, filename = temp_file, datatype = "FLT4S", overwrite = TRUE)
  return(raster_layer)
})

plot(gridMET_daily_vars[[1]][[1]], main = "Precipitation on 1998-01-01")
plot(st_geometry(shavers_prj), add = TRUE, border = "red", lwd = 2)

shavers_bounds <- vect(shavers_prj)


# ------------------------------------------------------------

# set CRS and clip a list of spatrasters to the watershed boundary
clip_to_watershed <- function(raster_list, bounds) {
  raster_list <- lapply(raster_list, function(raster_layer) {
    crs(raster_layer) <- "EPSG:4326"
    raster_layer
  })
  lapply(raster_list, function(raster_layer) {
    cropped <- crop(raster_layer, bounds)
    masked <- mask(cropped, bounds)
    return(masked)
  })
}

# loop through each variable and day to extract zonal stats
extract_daily_zonal_stats <- function(clipped_rasters, bounds) {
  all_stats <- data.frame()
  
  for (var_name in names(clipped_rasters)) {
    raster_layer <- clipped_rasters[[var_name]]
    
    for (i in 1:nlyr(raster_layer)) {
      date <- time(raster_layer)[i]
      
      cat("processing variable:", var_name, "on date:", as.character(date), "\n")
      
      extracted_values <- terra::extract(raster_layer[[i]], bounds)
      
      if (length(extracted_values) > 0 && !is.null(extracted_values[[1]])) {
        values <- unlist(extracted_values)
        stats <- data.frame(
          Year = year(date),
          Month = month(date),
          Day = day(date),
          Variable = var_name,
          Mean = mean(values, na.rm = TRUE),
          Median = median(values, na.rm = TRUE),
          Max = max(values, na.rm = TRUE),
          Min = min(values, na.rm = TRUE),
          Range = max(values, na.rm = TRUE) - min(values, na.rm = TRUE),
          SD = if(length(values) > 1) sd(values, na.rm = TRUE) else NA
        )
        all_stats <- dplyr::bind_rows(all_stats, stats)
      } else {
        cat("WARNING: no valid values for", var_name, "on date:", as.character(date), "\n")
      }
    }
  }
  return(all_stats)
}

# export zonal stats to individual CSVs per variable
export_stats_csv <- function(stats_df, output_dir) {
  for (var_name in unique(stats_df$Variable)) {
    var_data <- stats_df[stats_df$Variable == var_name, ]
    write.csv(var_data, file.path(output_dir, paste0("shavers_", var_name, ".csv")), row.names = FALSE)
  }
}


# --- gridMET: clip, extract zonal stats, export --------------------------------

shavers_daily_gridMET <- clip_to_watershed(gridMET_daily_vars, shavers_bounds)

plot(shavers_daily_gridMET[[1]][[1]], main = "Clipped Precipitation on 1998-01-01")
plot(shavers_bounds, add = TRUE, border = "red", lwd = 2)

shavers_gridMET_daily_stats <- extract_daily_zonal_stats(shavers_daily_gridMET, shavers_bounds)

for (var_name in unique(shavers_gridMET_daily_stats$Variable)) {
  var_data <- shavers_gridMET_daily_stats %>% filter(Variable == var_name)
  
  p <- ggplot(var_data, aes(x = as.Date(paste(Year, Month, Day, sep = "-")), y = Mean)) +
    geom_line(color = "blue") +
    labs(title = paste("Daily Mean of", var_name, "for Shavers Fork Watershed"),
         x = "Date", y = var_name) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
}

export_stats_csv(shavers_gridMET_daily_stats, "data/raw")


# --- MACA future climate data -------------------------------------------------

ESM2M <- "GFDL-ESM2M"
# MACA variables: huss, pr, rhsmax, rhsmin, rsds, tasmin, tasmax, vpd, uas, vas
MACA_vars <- c("huss", "pr", "rhsmax", "rhsmin", "rsds", "tasmin", "tasmax", "vpd", "uas", "vas")
rcp_scenarios = c("rcp85", "rcp45")
MACA_start_date <- "2024-01-01"
MACA_end_date <- "2099-12-31"

all_maca_vars_list <- list()

for (rcp in rcp_scenarios) {
  cat("starting extraction for scenario:", rcp, "\n")
  
  for (variable in MACA_vars) {
    cat("extracting:", variable, "(", rcp, ")\n")
    
    maca_data_extract <- getMACA(
      AOI = shavers_bounds,
      varname = variable,
      time = "day",
      model = ESM2M,
      scenario = rcp,
      startDate = MACA_start_date,
      endDate = MACA_end_date,
      verbose = TRUE,
      ID = NULL,
      dryrun = FALSE
    )
    
    if (is(maca_data_extract, "list") && length(maca_data_extract) == 1) {
      spat_raster <- maca_data_extract[[1]]
      
      if (is(spat_raster, "SpatRaster")) {
        list_name <- paste0(variable, "_", rcp)
        all_maca_vars_list[[list_name]] <- spat_raster
      } else {
        cat("Warning: expected SpatRaster but got", class(spat_raster), "for", variable, "\n")
      }
    } else {
      cat("Warning: unexpected extraction result for", variable, "\n")
    }
  }
}


# --- MACA: clip, extract zonal stats, export -----------------------------------

shavers_daily_maca_clipped <- clip_to_watershed(all_maca_vars_list, shavers_bounds)

plot(shavers_daily_maca_clipped[[2]][[1]], main = "Precipitation RCP 8.5 on 2024-01-01")
plot(shavers_bounds, add = TRUE, border = "red", lwd = 2)

shavers_MACA_daily_stats <- extract_daily_zonal_stats(shavers_daily_maca_clipped, shavers_bounds)

export_stats_csv(shavers_MACA_daily_stats, "data/raw")
