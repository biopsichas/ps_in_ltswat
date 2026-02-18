## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 1) Loading required libraries------------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Load multiple packages at once
packages <- c(
  "tidyverse", "lubridate", "data.table", "sf",
  "future.apply", "future", "DBI", "rlang",
  "rpostgis", "readxl", "RPostgres"
)

# Install missing packages automatically (optional)
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

# Load all packages
for(pkg in packages) library(pkg, character.only = TRUE)

# Source custom functions
source("function.R")

# SWATreadR from GitHub if not installed
if (!requireNamespace("SWATreadR", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("chrisschuerz/SWATreadR")
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## X) The rest------------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Connects to pg database and creates empty point source tables in point_sources_zero schema
## Defining SWAT LT database connection
lt_con <- tryCatch(
  dbConnect(
    RPostgres::Postgres(),
    dbname   = Sys.getenv("LT_DB_NAME"),
    host     = Sys.getenv("LT_DB_HOST"),
    port     = Sys.getenv("LT_DB_PORT"),
    user     = Sys.getenv("LT_DB_USER")
  ),
  error = function(e) stop("LT database connection failed: ", e$message)
)
## Check if successful
print(lt_con)
library(bit64)

## Data paths
gis_path <- "Data/GIS/"

## Loads existing small point source table
small_catch <- load_table(lt_con, "point_sources", "small_up2024")

## Checks the nutirent names. They should be ["SS","BOD7","NO3","NO2","NH4","NTOT","PO4","PTOT"]
unique(small_catch$nutrient)
names(small_catch)

## Updates
small_catch[small_catch$nutrient == "Ntot","nutrient"] <- "NTOT"
small_catch[small_catch$nutrient == "Ptot","nutrient"] <- "PTOT"

## Recheck
unique(small_catch$nutrient)

## Add cach_id
## Read basin shapefile
basins <- st_read(dsn = paste0(gis_path,"basin.gdb"), layer = "basin",
                  quiet = T) %>%
  select(c("cach_id"))

small_catch_up <- small_catch |>
  select(id, year, name, type, nutrient, concentration, unit, discharge, x, y) %>%
  st_as_sf(coords = c("x", "y"), crs = 3346, remove = FALSE) |>
  st_join(basins, left = TRUE) |>
  mutate(type = as.integer(type)) |>
  rename(catchmentid = cach_id,
         shape = geometry)|>
  select(id, year, name, type, nutrient, concentration, unit, discharge, x, y, shape, catchmentid)

## Write the table
table_constraint <- "
  id BIGINT,
  year INTEGER,
  name TEXT,
  type INTEGER,
  nutrient TEXT,
  concentration DOUBLE PRECISION,
  unit TEXT,
  discharge DOUBLE PRECISION,
  y DOUBLE PRECISION,
  x DOUBLE PRECISION,
  shape geometry(Point, 3346),
  catchmentid BIGINT
"
write_in_table(lt_con, "point_sources", "small_up2024_updated", small_catch_up, table_constraint = table_constraint)
