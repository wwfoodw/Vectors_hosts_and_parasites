# map

# Load libraries
library(ggplot2)
library(ggmap)
library(sf)


map_data <- read.csv("rangitahua_map_data.csv")
View(map_data)
map_data

library(basemaps)

# Bounding box with buffer
buffer <- 0.001
xmin <- min(map_data$Longitude) - buffer
xmax <- max(map_data$Longitude) + buffer
buffer <- 0.001
ymin <- min(map_data$Latitude) - buffer
ymax <- max(map_data$Latitude) + buffer

# Make the polygon geometry
polygon_geom <- st_sfc(
  st_polygon(list(matrix(c(
    xmin, ymin,
    xmax, ymin,
    xmax, ymax,
    xmin, ymax,
    xmin, ymin  # close the polygon
  ), ncol = 2, byrow = TRUE))),
  crs = 4326
)


# Wrap into an sf object called ext
ext <- st_sf(
  feature_type = "rectangle",
  geometry = polygon_geom
)


# view all available maps
get_maptypes()

set_defaults(map_service = "esri", map_type = "world_topo_map")
bm <- basemap_raster(ext, resolution = 2)

# load and return basemap map as class of choice, e.g. as image using magick:
basemap_magick(ext)

# or as plot:
basemap_plot(ext)

# or as ggplot2:
basemap_ggplot(ext)

## I have a map - now:

library(sf)
library(ggplot2)

# 1. Convert to sf
map_data_sf <- st_as_sf(map_data, coords = c("Longitude", "Latitude"), crs = 4326)

# 2. Transform to match basemap CRS (Web Mercator)
map_data_sf <- st_transform(map_data_sf, crs = 3857)

# 3. Plot
basemap_ggplot(ext) +
  geom_sf(data = map_data_sf[1,], color = "red", size = 3) +
  geom_sf(data = map_data_sf[c(2:6),], color = "green", size = 3) +
  geom_text(data = map_data_sf,
            aes(label = ID, geometry = geometry),
            stat = "sf_coordinates",
            size = 3, hjust = -0.1, vjust = -0.5) +
  theme(legend.position = "none")

###############
#
# Wider context
#
###############

library(sf)
library(ggplot2)
library(rnaturalearth)

# New Zealand outline
nz <- ne_countries(
  country = "New Zealand",
  scale = "medium",
  returnclass = "sf"
)

# Shift longitudes from -180–180 to 0–360
nz <- st_shift_longitude(nz)

# Bounding box around your Rangitāhua observations
buffer <- 0.1

rangitahua_box <- data.frame(
  xmin = min(map_data$Longitude) + 360 - buffer,
  xmax = max(map_data$Longitude) + 360 + buffer,
  ymin = min(map_data$Latitude) - buffer,
  ymax = max(map_data$Latitude) + buffer
)

ggplot() +
  geom_sf(
    data = nz,
    fill = "grey90",
    colour = "grey30"
  ) +
  geom_rect(
    data = rangitahua_box,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 181.8,
    y = -29.25,
    label = "Rangitāhua",
    colour = "red",
    hjust = 1
  ) +
  coord_sf(
    xlim = c(165, 183),
    ylim = c(-48, -17),
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  scale_x_continuous(
    breaks = c(165, 170, 175, 180),
    labels = function(x) {
      ifelse(
        x < 180,
        paste0(x, "°E"),
        ifelse(x == 180, "180°", paste0(360 - x, "°W"))
      ) 
    }
  ) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addRectangles(
    lng1 = xmin_r,
    lng2 = xmax_r,
    lat1 = ymin_r,
    lat2 = ymax_r,
    color = "red",
    weight = 3,
    fill = FALSE
  ) |>
  fitBounds(
    lng1 = 165,
    lat1 = -48,
    lng2 = 183,
    lat2 = -27
  )
