# Script to create a ggplot map in R

library(sf)
library(ggplot2)
library(ggspatial)

folder_data = "." # Folder of your shapefile
folder_results = "." # Folder where you want your output to be written

# Control the resolution and size of your figure
dpi_fig = 500
scale_fig = 0.8

# Link to a shapefile. You can create them in Google Earth
# See https://gis.stackexchange.com/questions/313662/how-can-i-generate-shapefile-from-google-earth
filename = "ErkenShp.kml"

# Data in decimal degrees (usually standard projection WGS84, crs = 4326)
se_water = st_read(file.path(folder_data, filename))

ggplot(se_water) +
  geom_sf()

### Add details and make figure nicer
# Type of terrain (just an example, but depending on your shape file, you may
# have different land types)
se_water$Type = c("water", rep("land", nrow(se_water) - 1L))

# Specific locations that you want to show on your map (again in decimal degrees)
data_points = data.frame(Latitude = c(59.83930583920924, 59.84297, 59.85513358640419, 59.85396582286886),
                         Longitude = c(18.62931534400336, 18.6354, 18.478677888268034, 18.650574136991295),
                         Location = c("Weather station & Nutrient sampling", "Thermistor chain", "Main inflow", "Main outflow"),
                         Num = 1:4)
data_points = st_as_sf(data_points, coords = c("Longitude", "Latitude"), crs = 4326) # crs argument denotes projection

ggplot() +
  geom_sf(data = se_water[se_water$Type == "water",], colour = "black", fill = "lightgrey") +
  geom_sf(data = se_water[se_water$Type == "land",], colour = "black", fill = "white") +
  geom_sf_label(data = data_points, aes(label = Num)) +
  # Or use geom_sf() to plot dots
  labs(x = element_blank(), y = element_blank()) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.text = element_blank()) +
  annotation_scale() +
  annotation_north_arrow(which_north = "true", location = "lb",
                         pad_x = unit(0.5, "cm"), pad_y = unit(1.0, "cm"),
                         style = north_arrow_fancy_orienteering())

ggsave(file.path(folder_results, "Map.png"), scale = scale_fig, dpi = dpi_fig)
