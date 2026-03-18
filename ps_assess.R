##------------------------------------------------------------------------------
## 1) Loading required libraries
##------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(data.table)
library(sf)
library(future.apply)
library(future)
if (!requireNamespace("SWATreadR", quietly = TRUE)) {
  devtools::install_github("chrisschuerz/SWATreadR")
}
# Source custom functions
source("function.R")

##------------------------------------------------------------------------------
## 2) Setting parameters to select and data paths
##------------------------------------------------------------------------------

## Path to the folder where LT SWAT system setups are stored
setup_path <- "E:/LIFE_AAA/swat_lt/Projects/Setup_2020_coarse/Watersheds/"

sc_zero <- "point_zero"

## Data paths
gis_path <- "Data/GIS/"

##------------------------------------------------------------------------------
## 3) Reading data
##------------------------------------------------------------------------------

## Read basin shapefile
basins <- st_read(dsn = paste0(gis_path,"basin.gdb"), layer = "basin",
                  quiet = T) %>%
  select(c("cach_id","Subbasin","Setup_name", "GRIDCODE"))

## Read reaches shae file
segments <- st_read(paste0(gis_path, "segments_coarse.shp"), quiet = T) %>%
  st_drop_geometry

## Preparing Subbasin Setup_name table
sstable <- basins %>% st_drop_geometry %>% select(Subbasin, Setup_name) %>% unique

##------------------------------------------------------------------------------
## 4) Reading SWAT output files and calculating concentration in parallel
##------------------------------------------------------------------------------

plan(multisession, workers = 15)

# Build list of paths first
all_tasks <- lapply(seq_len(nrow(sstable)), function(i) {
  list(
    base  = paste0(setup_path, sstable$Subbasin[i], "/", sstable$Setup_name[i],
                   "/", sc_zero, "/channel_sd_day.txt"),
    subb  = sstable$Subbasin[i],
    setup = sstable$Setup_name[i]
  )
})

## Run in parallel:
results <- future_lapply(all_tasks, function(task) {
  if (!file.exists(task$base)) return(NULL)
  read_model(task$base, "base",   task$subb, task$setup)
})

df <- rbindlist(Filter(Negate(is.null), results), use.names = TRUE)
##Close parallel workers
plan(sequential)

##------------------------------------------------------------------------------
## 5) Reading SWAT reservoir files (needed for the retention calculation)
##------------------------------------------------------------------------------

# Pre-generate all file paths (Vectorized)
f_paths <- paste0(setup_path, sstable$Subbasin, "/", sstable$Setup_name,
                  "/", sc_zero, "/hydrology.res")

# Use lapply to read files that exist and combine results
res <- lapply(f_paths, function(path) {if (file.exists(path)) read_res(path)
  else NULL}) |> rbindlist()
# saveRDS(res, "test/res.rds")

plan(multisession, workers = 15)

# Build list of paths first
all_tasks <- lapply(seq_len(nrow(sstable)), function(i) {
  list(
    base  = paste0(setup_path, sstable$Subbasin[i], "/", sstable$Setup_name[i],
                   "/", sc_zero, "/reservoir_day.txt"),
    subb  = sstable$Subbasin[i],
    setup = sstable$Setup_name[i]
  )
})

## Run in parallel:
results <- future_lapply(all_tasks, function(task) {
  if (!file.exists(task$base)) return(NULL)
  read_res2(task$base, task$subb, task$setup)
})

df_res <- rbindlist(Filter(Negate(is.null), results), use.names = TRUE)
##Close parallel workers
plan(sequential)
# saveRDS(df_res, "test/res_assessment2.rds")
