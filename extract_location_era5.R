# Function to extract time series from an ERA5-Land netcdf for a specific location
# ERA5 reanalysis forcing can be downloaded on
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview
# Unit info can be found on the ERA5 website and in the netcdf file

# Author: Jorrit Mesman

library(ncdf4)
library(lubridate)

extract_location_era5 = function(ncdf, latitude, longitude){
  ### Read
  nc = nc_open(ncdf)
  
  all_names = names(nc$var)
  
  time = ncvar_get(nc, "time") # in hours since 1900-01-01 00:00:00
  lats =  ncvar_get(nc, "latitude")
  lons = ncvar_get(nc, "longitude")
  
  # Get variables for right lat/lon (which are rounded to 1 decimal)
  round_lat = round(latitude, 1)
  round_lon = round(longitude, 1)
  
  lat_index = which(abs(lats - round_lat) < 0.001)
  lon_index = which(abs(lons - round_lon) < 0.001)
  
  if(length(lat_index) == 0L | length(lon_index) == 0L){
    stop("Requested lat/lon not in the netcdf file!")
  }
  
  ### Compile the data
  df = data.frame(datetime = as.POSIXct("1900-01-01") + hours(time))
  
  for(i in all_names){
    values_tmp = ncvar_get(nc, i)
    df[[i]] = values_tmp[lon_index, lat_index, 1:dim(values_tmp)[3]]
  }
  
  # c_u10 = ncvar_get(nc, "u10") # 10 metre U wind component (m/s)
  # c_v10 = ncvar_get(nc, "v10") # 10 metre V wind component (m/s)
  # c_d2m = ncvar_get(nc, "d2m") # 2 metre dewpoint temperature (K)
  # c_t2m = ncvar_get(nc, "t2m") # 2 metre temperature (K)
  # c_sp = ncvar_get(nc, "sp") # Surface pressure (Pa)
  # c_ssrd = ncvar_get(nc, "ssrd") # Surface solar radiation downwards (J/m2)
  # c_strd = ncvar_get(nc, "strd") # Surface thermal radiation downwards (J/m2)
  # c_tp = ncvar_get(nc, "tp") # Total precipitation (m)
  
  nc_close(nc)
  
  rm(nc, values_tmp)
  
  # IMPORTANT NOTES for radiation and precipitation:
  # ERA5 says "For the reanalysis, the accumulation period is over the 1 hour ending at the validity date and time"
  # "For the ensemble members, ensemble mean and ensemble spread, the processing period is over the 3 hours ending at the validity date and time."
  # "To convert to watts per square metre (W m-2 ), the accumulated values should be divided by the accumulation period expressed in seconds"
  # Use information from https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
  # Hourly data:
  # If h = 1 UTC, W_m2 = J_m2-h / 3600
  # Otherwise, W_m2 = (J_m2-h - J_m2-h-1) / 3600
  # Similar for rainfall
  
  return(df)
}
