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

##------------------------------------------------------------------------------
## 2) Setting parameters to select and data paths
##------------------------------------------------------------------------------

## Path to the folder where LT SWAT system setups are stored
setup_path <- "E:/LIFE_AAA/swat_lt/Projects/Setup_2020_coarse/Watersheds/"

sc_base <- "basecase"
sc_change <- "point_zero"

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

read_model <- function(f_path, sc_name, subb, setup) {
  dt <- as.data.table(SWATreadR::read_swat(f_path))
  dt[, year_month := sprintf("%04d-%02d", yr, mon)]
  # Compute concentrations
  dt[, tn_conc := ifelse(flo_out == 0, NA_real_,
                         ((orgn_out + nh3_out + no3_out + no2_out)*1e6)/(flo_out*86400*1000))]
  dt[, tp_conc := ifelse(flo_out == 0, NA_real_,
                         ((sedp_out + solp_out)*1e6)/(flo_out*86400*1000))]
  # Add identifiers once
  dt[, `:=`(
    Subbasin = subb,
    Setup_name = setup,
    Scenario = sc_name
  )]
  # Keep only necessary columns
  dt <- dt[, .(unit, year_month, flo_out, tn_conc, tp_conc,
               Subbasin, Setup_name, Scenario)]
  # Ultra-fast aggregation
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE),
           by = .(unit, year_month, Subbasin, Setup_name, Scenario)]

  return(dt)
}

# Build list of paths first
all_tasks <- lapply(seq_len(nrow(sstable)), function(i) {
  subb  <- sstable[i, "Subbasin"]
  setup <- sstable[i, "Setup_name"]
  f_base   <- paste0(setup_path, subb, "/", setup, "/", sc_base, "/channel_sd_day.txt")
  f_change <- paste0(setup_path, subb, "/", setup, "/", sc_change, "/channel_sd_day.txt")
  list(base = f_base, change = f_change, subb = subb, setup = setup)
})

## Run in parallel:
results <- future_lapply(all_tasks, function(task) {
  if (!file.exists(task$base)) return(NULL)
  base_dt   <- read_model(task$base, "base",   task$subb, task$setup)
  change_dt <- read_model(task$change, "change", task$subb, task$setup)
  rbindlist(list(base_dt, change_dt))
})

df <- rbindlist(Filter(Negate(is.null), results), use.names = TRUE)
##Close parallel workers
plan(sequential)

##------------------------------------------------------------------------------
## 5) Cleaning up modeling data
##------------------------------------------------------------------------------

## Just to add cach_id info
basins_df <- basins %>%
  st_drop_geometry %>%
  select(cach_id, Subbasin, Setup_name, GRIDCODE)

# Convert to a named list for fast lookup by 'id'
segments_lookup <- segments %>%
  select(id, wbriver_co, wlake_wb)

## Joining SWAT output data with basin information, cleaning names, filtering
df_fix <- df %>%
  rename(GRIDCODE = unit) %>%
  left_join(basins_df, by = c("Subbasin", "Setup_name", "GRIDCODE")) %>%
  left_join(segments_lookup, by = c("cach_id" = "id")) %>%
  select(-c("GRIDCODE", "cach_id")) %>%
  mutate(tn_good = ifelse(tn_conc <= 3, TRUE, FALSE),
         tp_good = ifelse(tp_conc <= 0.14, TRUE, FALSE)) %>%
  mutate(month = as.integer(substr(year_month, 6, 7))) %>%
  select(-year_month) %>%
  group_by(Scenario, Subbasin, Setup_name, wbriver_co, wlake_wb, month) %>%
  summarise(flo_out = mean(flo_out, na.rm = TRUE),
            tn_conc = mean(tn_conc, na.rm = TRUE),
            tp_conc = mean(tp_conc, na.rm = TRUE),
            tn_good_perc = mean(tn_good, na.rm = TRUE) * 100,
            tp_good_perc = mean(tp_good, na.rm = TRUE) * 100, .groups = "drop") %>%
  filter(!is.na(wbriver_co) | !is.na(wlake_wb)) %>%
  select(-c("Subbasin", "Setup_name")) %>%
  setDT

dt_max <- df_fix[, .SD[which.max(flo_out)], by = .(Scenario, wbriver_co, wlake_wb, month)]

##------------------------------------------------------------------------------
## 6) Calculate indicators
##------------------------------------------------------------------------------

indicators <- dt_max %>%
  group_by(Scenario, wbriver_co, wlake_wb) %>%
  summarise(
    avg_tn_conc = mean(tn_conc, na.rm = TRUE),
    avg_tp_conc = mean(tp_conc, na.rm = TRUE),
    perc_tn_good = mean(tn_good_perc, na.rm = TRUE),
    perc_tp_good = mean(tp_good_perc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Scenario,
    values_from = c(avg_tn_conc, avg_tp_conc, perc_tn_good, perc_tp_good),
    names_sep = "_"
  ) %>%
  mutate(imp_tn = ifelse(avg_tn_conc_base > 3 & avg_tn_conc_change <= 3, TRUE, FALSE),
         imp_tp = ifelse(avg_tp_conc_base > 0.14 & avg_tp_conc_change <= 0.14, TRUE, FALSE))

important_wb <- indicators %>%
  filter(imp_tn == TRUE | imp_tp == TRUE)
