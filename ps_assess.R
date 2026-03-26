## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 1) Loading required libraries ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(tidyverse)
library(lubridate)
library(data.table)
library(sf)
library(future.apply)
library(future)
library(openxlsx)
library(rmapshaper)
if (!requireNamespace("SWATreadR", quietly = TRUE)) {
  devtools::install_github("chrisschuerz/SWATreadR")
}
# Source custom functions
source("function.R")

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 2) Setting parameters to select and data paths ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Read model if TRUE, otherwise read from RDS (saved results)
read_model_files <- FALSE

## Simplify geometry if TRUE for leaflet or read from RDS (saved results)
simplify_geometry <- TRUE

## Remove existing temp_folder if it exists (to ensure clean start)
clean_start <- FALSE

## Temp folder for saving model outputs
temp_folder <- "Temp"

## Path to the folder where LT SWAT system setups are stored
setup_path <- "E:/LIFE_AAA/swat_lt/Projects/Setup_2020_coarse/Watersheds/"

## Subdirectory names for the scenario to be used for the assessment
sc_zero <- "point_zero"

## Data paths
data_path <- "Data/"

## GIS data
gis_path <- paste0(data_path, "GIS/")
## Point source loads

## Water body to cach_id relationship, and problematic water bodies that need
## to be checked manually
wb_rel <- paste0(data_path, "WB_representative_SWAT_segments_final.xlsx")

## Point source loads
ps_data <- paste0(data_path, "PS/")

## Retention coefficients for rivers (per %/km)
river_N_retention <- 0.00024
river_P_retention <- 0.00045

## Defining number of cores for parallel processing (adjust based on your system)
num_cores <- 15

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 3) Reading data ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Manage the temp_folder directory
if(clean_start){
  if (dir.exists(temp_folder)) {
    unlink(temp_folder, recursive = TRUE) # Delete existing folder
    paste0("Existing folder '", temp_folder, "' deleted for a clean start.")
  }
  dir.create(temp_folder, recursive = TRUE) # Create fresh folder
} else {
  if (!dir.exists(temp_folder)) {
    dir.create(temp_folder, recursive = TRUE) # Create folder if it doesn't exist
    paste0("Folder '", temp_folder, "' created.")
  } else {
    paste0("Folder '", temp_folder, "' already exists. Some files may be overwritten.")
  }
}

## Read basin shapefile
basins_sf <- st_read(dsn = paste0(gis_path,"basin.gdb"), layer = "basin",
                     quiet = T) |>
  select(c("cach_id","Subbasin","Setup_name", "GRIDCODE"))

## Create a simplified table for easier joins later on
basins <- basins_sf |> st_drop_geometry() |>
  select(cach_id, Subbasin, Setup_name, GRIDCODE)

## Read reaches shae file
segments_sf <- st_read(paste0(gis_path, "segments_coarse.shp"), quiet = T)

if(simplify_geometry){
  # Check size in Megabytes basin
  print(paste0("Basin file is ", format(object.size(basins_sf), units = "Mb")))

  # Check size in Megabytes segments
  print(paste0("Basin file is ", format(object.size(segments_sf), units = "Mb")))

  basins_simp <- st_make_valid(basins_sf) |>
    st_simplify(preserveTopology = TRUE, dTolerance = 10) |>
    ms_simplify(
      keep = 0.05,          # 5% of original detail
      snap = TRUE,          # Fuses boundaries together
      snap_interval = 0.001, # Adjust this if gaps persist
      keep_shapes = TRUE
    ) |> st_transform(4326)

  segments_simp <- st_simplify(segments_sf, preserveTopology = TRUE, dTolerance = 10) |>
    st_transform(4326)

  # mapview::mapview(segments_simp) + mapview::mapview(basins_simp)

  print(paste0("Basin file is ", format(object.size(basins_simp), units = "Mb")))
  print(paste0("Basin file is ", format(object.size(basins_simp), units = "Mb")))

  saveRDS(basins_simp, file = file.path(temp_folder, "basins_simp.rds"))
  saveRDS(segments_simp, file = file.path(temp_folder, "segments_simp.rds"))
} else {
  ## Reading rds files if they exist,
  ## otherwise will need to be created by setting simplify_geometry = TRUE
  basins_simp <- readRDS(file.path(temp_folder, "basins_simp.rds"))
  segments_simp <- readRDS(file.path(temp_folder, "segments_simp.rds"))
  print(paste0("Simplified basin and segment files loaded from RDS"))
}


## Preparing Subbasin Setup_name table
sstable <- basins_sf %>% st_drop_geometry %>% select(Subbasin, Setup_name) %>%
  unique

## Reading tables
## Transfers table: This table contains information about water transfers between
## catchments, which can affect the flow and pollutant loads in the river network.
transfer <- read.csv(paste0(data_path, "transfers.csv"),
                     header = FALSE,
                     fileEncoding = "UTF-8") [,-1]|>
  setNames(c("id", "name", "outletid", "flowto", "multiplier")) |>
  select(name, outletid, flowto, multiplier)

## Catchment coarse info: This table contains information about the catchments,
## such as their type, area, and flow relationships. This information is crucial
## for understanding the characteristics of each catchment and how they contribute
## to the overall river network (catchments_coarse layer attributes).

catch_coarse_info <- read.csv("Data/catchment_coarse.csv", header = F) |>
  setNames(c("id", "type", "segmentid", "lakegid", "kadastroid", "kadastroid_lake",
             "flow_to", "segmentto", "outletid", "area", "addedarea", "inflowarea"))

## Reading the relationship between water bodies and catchment IDs, as well as a list
## of problematic water bodies that need to be checked manually.
wb_to_cach_id <- read.xlsx(wb_rel, sheet = "cach_id_to_wb")
wb_problematic <- read.xlsx(wb_rel, sheet = "problem_wb")

## Reading point source loads: This dataset contains information about
## the pollutant loads from point sources (e.g., wastewater treatment plants,
## industrial discharges) into the river network. Reading is done into the list
## form the ps_data folder, where each file corresponds
## to a different point source scenario.

## List all files in the point source data folder
ps_files <- list.files(ps_data, pattern = "\\.xlsx$", full.names = TRUE)
ps_names <- tools::file_path_sans_ext(basename(ps_files))
ps_lst <- list()
for(i in seq_along(ps_files)) {
  ps_lst[[ps_names[i]]] <- read.xlsx(ps_files[i], sheet = 1, check.names = FALSE)[c(1:7)] |>
    setNames(c("Išleistuvo kodas", "Pavadinimas", "Nuotėkų kiekis 1000 m3/metus",
               "Bendrasis azotas (kg/metus)", "Bendrasis fosforas (kg/metus)", "X", "Y")) |>
    mutate(across(c(1, 6, 7), as.integer),
           across(c(3, 4, 5), as.numeric),
           across(c(3, 4, 5), ~round(.,1))) |>
    filter(X > 0 & Y > 0) |> # Filter out rows with invalid coordinates
    st_as_sf(coords = c("X", "Y"), crs = 3346, remove = FALSE) |>
    st_join(basins_sf["cach_id"], left = TRUE) |>
    st_drop_geometry() |>
    select(cach_id, everything())
  print(paste0("Point source data in ", ps_names[i], " was read."))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4) Reading SWAT output files and calculating concentration in parallel ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if(read_model_files){
  ## Define parallel plan
  # Note: Ensure 'workers' doesn't exceed your physical cores to avoid overhead
  plan(multisession, workers = num_cores)

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
    # Assuming read_model is a custom function you've defined elsewhere
    read_model(task$base, "base", task$subb, task$setup)
  })

  ## Combine results, filtering out any NULLs
  df_mod <- rbindlist(Filter(Negate(is.null), results), use.names = TRUE)

  ## Close parallel workers
  plan(sequential)

  ## --- Save to temp_folder .rds file ---
  # Save the result
  saveRDS(df_mod, file = file.path(temp_folder, "df_mod.rds"))

  message("Model files processed and saved to: ", temp_folder)
} else {
  ## If not reading raw files, load the previously saved .rds data
  processed_file <- file.path(temp_folder, "df_mod.rds")

  if (file.exists(processed_file)) {
    message("Loading previously processed model data from: ", processed_file)
    df_mod <- readRDS(processed_file)
  } else {
    stop("The processed file 'df_mod.rds' was not found in ", temp_folder,
         ". Please set read_model_files = TRUE to generate it.")
  }
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 5) Reading SWAT reservoir files (needed for the retention calculation) ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if(read_model_files){
  ## Define parallel plan
  plan(multisession, workers = num_cores)

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

  ## Combine results, filtering out any NULLs
  df_res <- rbindlist(Filter(Negate(is.null), results), use.names = TRUE)

  ## Close parallel workers
  plan(sequential)

  ## --- Save to temp_folder ---
  # Ensure the directory exists (keeping existing folder if created by previous block)
  if (!dir.exists(temp_folder)) {
    dir.create(temp_folder, recursive = TRUE)
  }

  saveRDS(df_res, file = file.path(temp_folder, "df_res.rds"))
  message("Reservoir files processed and saved.")

} else {
  ## Load the previously saved reservoir data
  res_file <- file.path(temp_folder, "df_res.rds")

  if (file.exists(res_file)) {
    df_res <- readRDS(res_file)
    message("Loaded cached reservoir data from: ", res_file)
  } else {
    stop("The file 'df_res.rds' is missing. Please set read_model_files = TRUE.")
  }
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 6) Cleaning all the information ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Cleaning and summarizing model output data
df_mod_cl <- df_mod |>
  left_join(basins,
            by = c("unit" = "GRIDCODE", "Subbasin", "Setup_name")) |>
  select(cach_id, flo_out, tn_conc, tp_conc) |>
  group_by(cach_id) |>
  summarise(across(c(flo_out, tn_conc, tp_conc), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") |>
  # Ensure all IDs are present
  right_join(basins |> select(cach_id), by = "cach_id") |>
  # Efficiently replace NAs
  mutate(across(c(flo_out, tn_conc, tp_conc), ~ replace_na(.x, 0)))

## Cleaning and summarizing reservoir data
df_res_cl <- df_res |>
  left_join(basins, by = c("unit" = "GRIDCODE", "Subbasin", "Setup_name")) |>
  select(cach_id, flo_stor, flo_out, tn_in, tn_out, tp_in, tp_out) |>
  group_by(cach_id) |>
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  # Calculate metrics on the averaged values
  mutate(
    wrt_days = if_else(flo_out > 0, flo_stor / flo_out, 0),
    l_retN   = pmax(0, 1 - (tn_out / tn_in)), # pmax is faster/cleaner for clamping to 0
    l_retP   = 1 - (tp_out / tp_in)
  ) |>
  # Final alignment with the full basin list
  right_join(basins |> select(cach_id), by = "cach_id") |>
  select(cach_id, wrt_days, l_retN, l_retP) |>
  # Clean up NAs
  mutate(across(c(wrt_days, l_retN, l_retP), ~ replace_na(.x, 0)))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 7) MAIN PART: Building the river map with all the information. ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 8) Preparing the data for presentation ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

