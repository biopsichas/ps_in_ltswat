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

## Retention coefficients for rivers (per %/km)
river_N_retention <- 0.00024
river_P_retention <- 0.00045

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

## Reading tables
transfer <- read.csv("Data/transfers.csv", header = F) [,-1]|>
  setNames(c("id", "name", "outletid", "flowto", "multiplier")) |>
  select(name, outletid, flowto, multiplier)

catch_coarse_info <- read.csv("Data/catchment_coarse.csv", header = F) |>
  setNames(c("id", "type", "segmentid", "lakegid", "kadastroid", "kadastroid_lake",
             "flow_to", "segmentto", "outletid", "area", "addedarea", "inflowarea"))

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

df_mod <- rbindlist(Filter(Negate(is.null), results), use.names = TRUE)
##Close parallel workers
plan(sequential)

##------------------------------------------------------------------------------
## 5) Reading SWAT reservoir files (needed for the retention calculation)
##------------------------------------------------------------------------------

# ## Saving the information from reservoir files
# # Pre-generate all file paths (Vectorized)
# f_paths <- paste0(setup_path, sstable$Subbasin, "/", sstable$Setup_name,
#                   "/", sc_zero, "/hydrology.res")
#
# # Use lapply to read files that exist and combine results
# res <- lapply(f_paths, function(path) {if (file.exists(path)) read_res(path)
#   else NULL}) |> rbindlist()
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
# saveRDS(df_res, "test/res_assessment2.rds")
##Close parallel workers
plan(sequential)

##------------------------------------------------------------------------------
## 6) Cleaning all the information
##------------------------------------------------------------------------------

df_mod <- df_mod |>
  left_join(basins |> st_drop_geometry(), by =c("unit" = "GRIDCODE", "Subbasin", "Setup_name")) |>
  select(cach_id, flo_out, tn_conc, tp_conc) |>
  group_by(cach_id) |>
  summarise_all(mean, na.rm = TRUE) |>
  right_join(basins[c("cach_id")] |> st_drop_geometry(), by = "cach_id") |>
  replace_na(list(flo_out = 0, tn_conc = 0, tp_conc = 0))

df_res <- df_res |>
  left_join(basins |> st_drop_geometry(), by =c("unit" = "GRIDCODE", "Subbasin", "Setup_name")) |>
  select(cach_id, year_month, flo_stor, flo_out, tn_in, tn_out, tp_in, tp_out) |>
  select(-year_month) |>
  group_by(cach_id) |>
  summarise_all(mean, na.rm = TRUE) |>
  mutate(wrt_days = ifelse(flo_out > 0, flo_stor/ flo_out, 0),#Calculate WRT in days based on average volume|>
         l_retN = 1 - (tn_out/tn_in),
         l_retP = 1 - (tp_out/tp_in),
         l_retN = ifelse(l_retN < 0, 0, l_retN)) |> # Set negative retention to 0
  right_join(basins[c("cach_id")] |> st_drop_geometry(), by = "cach_id") |>
  select(cach_id, wrt_days, l_retN, l_retP) |>
  replace_na(list(wrt_days = 0, l_retN = 0, l_retP = 0))

df_ps <- read.csv("point_source_loads.csv",
                  check.names = FALSE,
                  fileEncoding = "UTF-8")|>
  select(c("Išleistuvo kodas", "Teršalo pavadinimas", "Tiesioginis VT cach_id",  "Teršalo kiekis išleidžiamose nuotekose (kg/metus)")) |>
  filter(`Teršalo pavadinimas` %in% c("Bendrasis azotas", "Bendrasis fosforas")) |>
  setNames(c("ps_code", "populant", "cach_id", "load_kg_per_year")) |>
  mutate(cach_id = as.integer(cach_id),
         populant = ifelse(populant == "Bendrasis azotas", "TN", "TP"),
         load_kg_per_year = as.numeric(load_kg_per_year)) |>
  pivot_wider(
    id_cols = c(cach_id, ps_code), # This tells R: "one row per code per catchment"
    names_from = populant,
    values_from = load_kg_per_year,
    values_fn = sum
  ) |>
  replace_na(list(TN = 0, TP = 0))

## Converting GIS files for the leaflet
basins_4326 <- st_transform(basins, 4326)
segments_4326 <- st_read(paste0(gis_path, "segments_coarse.shp"), quiet = T) |>
  st_transform(4326)
ps_info <- readRDS("Data/ps_sf_info.rds")

##------------------------------------------------------------------------------
## 7) MAIN PART: Building the river map with all the information, starting from upstream and moving downstream
##------------------------------------------------------------------------------

# 1. Prepare info table
basins_info <- basins |>
  st_drop_geometry() |>
  left_join(catch_coarse_info, by = c("cach_id" = "id")) |>
  select(cach_id, type, segmentid, flow_to, segmentto, outletid, area, addedarea, inflowarea) |>
  left_join(segments[,c("id", "length")], by = c("segmentid" = "id")) |>
  mutate(len = ifelse(is.na(length), 0, length/1000),
         r_retN = 1-(len * river_N_retention),
         r_retP = 1-(len * river_P_retention))

# 2. Identify initial upstream IDs (IDs that don't receive flow from anywhere else)
upstreams_ids <- setdiff(basins_info$cach_id, basins_info$flow_to)

# 3. Filter out IDs that exist in the transfer table's flowto column
# This replaces your 'x' logic entirely in one clean step
upstreams_ids <- upstreams_ids[!(upstreams_ids %in% transfer$flowto)]


river_map <- lapply(upstreams_ids, function(id) {
  ps <- df_ps[df_ps$cach_id == id,]
  df <- df_mod[df_mod$cach_id == id,]
  # Return a nested list: one part for the area, one part for the inflows
  inflow = data.frame(
    cach_id = id,
    fraction = 1,
    area = basins_info$area[match(id, basins_info$cach_id)],
    added_area = basins_info$addedarea[match(id, basins_info$cach_id)],
    inflow_area = basins_info$inflowarea[match(id, basins_info$cach_id)],
    len = basins_info$length[match(id, basins_info$cach_id)],
    r_retN = basins_info$r_retN[match(id, basins_info$cach_id)],
    r_retP = basins_info$r_retP[match(id, basins_info$cach_id)],
    l_retN = 1 - df_res$l_retN[match(id, df_res$cach_id)],
    l_retP = 1 - df_res$l_retP[match(id, df_res$cach_id)],
    zero_N_load = df$flo_out * df$tn_conc * 31536,
    zero_P_load = df$flo_out * df$tp_conc * 31536,
    stringsAsFactors = FALSE)
  if(nrow(ps)>0){
    ps_load = ps |>
      mutate(TN = TN * inflow$r_retN * inflow$l_retN,
             TP = TP * inflow$r_retP * inflow$l_retP,
             TN_conc_added = TN / (df$flo_out * 31536),
             TP_conc_added = TP / (df$flo_out * 31536))
  } else {
    ps_load = data.frame(
      cach_id = integer(),
      ps_code = integer(),
      TN = numeric(),
      TP = numeric(),
      TN_conc_added = numeric(),
      TP_conc_added = numeric(),
      stringsAsFactors = FALSE
    )
  }

  list(
    inflow = inflow,
    ps_load = ps_load
  )
})

names(river_map) <- upstreams_ids

# Get a list of all unique basin IDs that need to be processed
to_do_ids <- setdiff(unique(basins_info$cach_id), names(river_map))

# Keep looping while there are IDs in basins_info not yet in river_map
while(length(to_do_ids) > 0) {

  # Check how many we have at the start of this pass (to prevent infinite loops)
  starting_count <- length(to_do_ids)

  for(id in to_do_ids) {

    # Get all potential upstream contributors for this 'id'
    inflow_ids <- basins_info$cach_id[which(basins_info$flow_to == id)]
    inflow_tranfer_id <- NULL

    if(id %in% transfer$flowto) {
      outlet <- transfer$outletid[match(id, transfer$flowto)]
      inflow_tranfer_id <- basins_info$cach_id[which(basins_info$outletid == outlet)]
      inflow_ids <- c(inflow_ids, inflow_tranfer_id)
      multiplier <- transfer$multiplier[match(id, transfer$flowto)]
    }
    # CRITICAL: Only proceed if ALL upstream contributors are now in river_map
    if(all(as.character(inflow_ids) %in% names(river_map))) {
      ps <- df_ps[df_ps$cach_id == id,]
      df <- df_mod[df_mod$cach_id == id,]

      inflow <- data.frame(
        cach_id = id,
        fraction = 1,
        area = basins_info$area[match(id, basins_info$cach_id)],
        added_area = basins_info$addedarea[match(id, basins_info$cach_id)],
        inflow_area = basins_info$inflowarea[match(id, basins_info$cach_id)],
        len = basins_info$length[match(id, basins_info$cach_id)],
        r_retN = basins_info$r_retN[match(id, basins_info$cach_id)],
        r_retP = basins_info$r_retP[match(id, basins_info$cach_id)],
        l_retN = 1 - df_res$l_retN[match(id, df_res$cach_id)],
        l_retP = 1 - df_res$l_retP[match(id, df_res$cach_id)],
        zero_N_load = df$flo_out * df$tn_conc * 31536,
        zero_P_load = df$flo_out * df$tp_conc * 31536,
        stringsAsFactors = FALSE
      )
      inflow_br <- inflow[c("r_retN", "r_retP", "l_retN", "l_retP")]

      if(nrow(ps)>0){
        ps_load = ps |>
          mutate(TN = TN * inflow$r_retN * inflow$l_retN,
                 TP = TP * inflow$r_retP * inflow$l_retP,
                 TN_conc_added = TN / (df$flo_out * 31536),
                 TP_conc_added = TP / (df$flo_out * 31536))
      } else {
        ps_load = data.frame(
          cach_id = integer(),
          ps_code = integer(),
          TN = numeric(),
          TP = numeric(),
          TN_conc_added = numeric(),
          TP_conc_added = numeric(),
          stringsAsFactors = FALSE
        )
      }
      for(in_id in inflow_ids) {
        if(!is.null(inflow_tranfer_id) && in_id %in% inflow_tranfer_id){
          tmp <- river_map[[as.character(in_id)]]$inflow
          tmp$fraction <- tmp$fraction * multiplier

          if(nrow(river_map[[as.character(in_id)]]$ps_load) > 0){
            ps_load_tmp <- river_map[[as.character(in_id)]]$ps_load |>
              mutate(TN = TN * multiplier * inflow_br $r_retN * inflow_br$l_retN,
                     TP = TP * multiplier * inflow_br $r_retP * inflow_br$l_retP,
                     TN_conc_added = TN / (df$flo_out * 31536),
                     TP_conc_added = TP / (df$flo_out * 31536))
          }
        } else {
          tmp <- river_map[[as.character(in_id)]]$inflow

          if(nrow(river_map[[as.character(in_id)]]$ps_load) > 0){
            ps_load_tmp <- river_map[[as.character(in_id)]]$ps_load |>
              mutate(TN = TN * inflow_br$r_retN * inflow_br$l_retN,
                     TP = TP * inflow_br$r_retP * inflow_br$l_retP,
                     TN_conc_added = TN / (df$flo_out * 31536),
                     TP_conc_added = TP / (df$flo_out * 31536))
          }
        }

        inflow <- rbind(inflow, tmp)
        if(exists("ps_load_tmp")){
          ps_load <- rbind(ps_load, ps_load_tmp)
          rm(ps_load_tmp)
        }
      }
      river_map[[as.character(id)]]$inflow <- inflow
      river_map[[as.character(id)]]$ps_load <- ps_load
      to_do_ids <- setdiff(to_do_ids, id)
    }
  }

  # Safety Break: If a whole loop finishes and we didn't add anything new,
  # something is wrong with the network logic (e.g., a circular reference)
  if(length(to_do_ids) == starting_count) {
    warning("Loop stopped: No new basins could be processed. Check for gaps in your network.")
    break
  }
}

# saveRDS(river_map, "test/river_map.rds")

##------------------------------------------------------------------------------
## 8) Preparing the data for presentation (rounding, calculating shares, etc.)
##------------------------------------------------------------------------------

for(id in names(river_map)){
  river_map[[id]]$inflow <- river_map[[id]]$inflow |>
    mutate(area = round(area/1000000, 1),
           added_area = round(added_area/1000000, 1),
           inflow_area = round(inflow_area/1000000, 1),
           len = round(len/1000, 1),
           r_retN = round(r_retN, 3),
           r_retP = round(r_retP, 3),
           l_retN = round(l_retN, 3),
           l_retP = round(l_retP, 3),
           zero_N_load = round(zero_N_load/1000, 1),
           zero_P_load = round(zero_P_load/1000, 1))
  ps_n_sum <- river_map[[id]]$ps_load$TN |> sum()
  ps_p_sum <- river_map[[id]]$ps_load$TP |> sum()
  river_map[[id]]$ps_load <- river_map[[id]]$ps_load |>
    mutate(TN_proc = round(100*(TN/ps_n_sum), 1),
           TP_proc = round(100*(TP/ps_p_sum), 1),
           TN = round(TN/1000, 3),
           TP = round(TP/1000, 3),
           TN_conc_added = round(TN_conc_added, 3),
           TP_conc_added = round(TP_conc_added, 3))
}

