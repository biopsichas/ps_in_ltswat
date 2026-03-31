## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 1) Loading required libraries ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(tidyverse)
library(lubridate)
library(data.table)
library(sf)
library(future)
library(future.apply)
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

## Simplify geometry if TRUE for leaflet or FALSE to read from RDS (saved results)
simplify_geometry <- FALSE

## Remove existing temp_folder if it exists (to ensure clean start)
clean_start <- FALSE

## Save calculated river_map as RDS file in temp_folder (TRUE) or not (FALSE)
save_river_map <- TRUE

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
ps_data <- paste0(data_path, "PS/")

## Retention coefficients for rivers (per %/km)
river_N_retention <- 0.00024
river_P_retention <- 0.00045

## Defining number of cores for parallel processing (adjust based on your system)
## Only used if read_model_files = TRUE,
## otherwise will read from RDS which is much faster.
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
segments <- segments_sf |> st_drop_geometry()

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

# 4. Initialize river_map with the upstream IDs, including their inflow and point source loads
river_map <- lapply(upstreams_ids, function(id) {
  # Basic catchment data
  df <- df_mod_cl[df_mod_cl$cach_id == id,]
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
    l_retN = 1 - df_res_cl$l_retN[match(id, df_res_cl$cach_id)],
    l_retP = 1 - df_res_cl$l_retP[match(id, df_res_cl$cach_id)],
    zero_N_load = df$flo_out * df$tn_conc * 31536,
    zero_P_load = df$flo_out * df$tp_conc * 31536,
    stringsAsFactors = FALSE)
  # 2. Iterate through each dataframe in ps_lst
  # This creates a list of ps_load dataframes for this specific ID
    ps_load_list <- lapply(names(ps_lst), function(nm) {
      ps_subset <- ps_lst[[nm]][ps_lst[[nm]]$cach_id == id & !is.na(ps_lst[[nm]]$cach_id), ]
      if(nrow(ps_subset) > 0) {
        # Apply your logic to the specific scenario
        ps_res <- ps_subset |>
          mutate(
            TN = `Bendrasis azotas (kg/metus)` * inflow$r_retN * inflow$l_retN,
            TP = `Bendrasis fosforas (kg/metus)` * inflow$r_retP * inflow$l_retP,
            TN_conc_added = TN / (df$flo_out * 31536),
            TP_conc_added = TP / (df$flo_out * 31536)
          )
      } else {
        # Empty dataframe with correct columns if no match found for this ID
        ps_res <- data.frame(
          cach_id = integer(),
          TN = numeric(),
          TP = numeric(),
          TN_conc_added = numeric(),
          TP_conc_added = numeric(),
          stringsAsFactors = FALSE
        )
      }
      return(ps_res)
    })

  # Name the internal list elements so you can call river_map[[1]]$ps_load$scenario_name
  names(ps_load_list) <- names(ps_lst)

  list(
    inflow = inflow,
    ps_load = ps_load_list
  )
})

## Assign names to river_map based on the upstream IDs for easy access
names(river_map) <- upstreams_ids

## Check the size of the river_map object in memory
print(paste0("Upstreams calculated 'river_map' object size ",
             format(object.size(river_map), units = "MB")))

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
      df <- df_mod_cl[df_mod_cl$cach_id == id,]

      # 1. Standard Inflow Calculation (remains a single dataframe)
      inflow <- data.frame(
        cach_id = id,
        fraction = 1,
        area = basins_info$area[match(id, basins_info$cach_id)],
        added_area = basins_info$addedarea[match(id, basins_info$cach_id)],
        inflow_area = basins_info$inflowarea[match(id, basins_info$cach_id)],
        len = basins_info$length[match(id, basins_info$cach_id)],
        r_retN = basins_info$r_retN[match(id, basins_info$cach_id)],
        r_retP = basins_info$r_retP[match(id, basins_info$cach_id)],
        l_retN = 1 - df_res_cl$l_retN[match(id, df_res_cl$cach_id)],
        l_retP = 1 - df_res_cl$l_retP[match(id, df_res_cl$cach_id)],
        zero_N_load = df$flo_out * df$tn_conc * 31536,
        zero_P_load = df$flo_out * df$tp_conc * 31536,
        stringsAsFactors = FALSE
      )
      inflow_br <- inflow[c("r_retN", "r_retP", "l_retN", "l_retP")]

      # 2. Iterate through each point source scenario, the current cach_id
      ps_load_list <- lapply(names(ps_lst), function(nm) {
        ps_subset <- ps_lst[[nm]][ps_lst[[nm]]$cach_id == id & !is.na(ps_lst[[nm]]$cach_id), ]
        if(nrow(ps_subset) > 0) {
          ps_load_acc <- ps_subset |>
            mutate(TN = `Bendrasis azotas (kg/metus)` * inflow$r_retN * inflow$l_retN,
                   TP = `Bendrasis fosforas (kg/metus)` * inflow$r_retP * inflow$l_retP,
                   TN_conc_added = TN / (df$flo_out[1] * 31536),
                   TP_conc_added = TP / (df$flo_out[1] * 31536))
        } else {
          ps_load_acc <- data.frame(
            "cach_id" = numeric(),
            "Išleistuvo kodas" = integer(),
            "Pavadinimas" = character(),
            "Nuotėkų kiekis 1000 m3/metus" = numeric(),
            "Bendrasis azotas (kg/metus)" = numeric(),
            "Bendrasis fosforas (kg/metus)" = numeric(),
            "X" = integer(),
            "Y" = integer(),
            "TN" = numeric(),
            "TP" = numeric(),
            "TN_conc_added" = numeric(),
            "TP_conc_added" = numeric(),
            stringsAsFactors = FALSE,
            check.names = FALSE # CRITICAL: prevents R from changing spaces to dots
          )
        }
        # Pull and aggregate loads from all upstream IDs for THIS scenario
        for(in_id in inflow_ids) {
          up_ps <- river_map[[as.character(in_id)]]$ps_load[[nm]]
          if(nrow(up_ps) > 0) {
            # Apply transfer multiplier if applicable
            m <- if(!is.null(inflow_tranfer_id) && in_id %in% inflow_tranfer_id) multiplier else 1

            ## Apply retention and concentration calculations to the upstream point source loads
            ps_up_processed <- up_ps |>
              mutate(TN = TN * m * inflow$r_retN * inflow$l_retN,
                     TP = TP * m * inflow$r_retP * inflow$l_retP,
                     TN_conc_added = TN / (df$flo_out[1] * 31536),
                     TP_conc_added = TP / (df$flo_out[1] * 31536))

            ps_load_acc <- rbind(ps_load_acc, ps_up_processed)
          }
        }
        return(ps_load_acc)
      })
      names(ps_load_list) <- names(ps_lst)

      # 3. Handle physical inflow geometry (Area accumulation)
      for(in_id in inflow_ids) {
        tmp <- river_map[[as.character(in_id)]]$inflow
        if(!is.null(inflow_tranfer_id) && in_id %in% inflow_tranfer_id) tmp$fraction <- tmp$fraction * multiplier
        inflow <- rbind(inflow, tmp)
      }

      # Save results to the map
      river_map[[as.character(id)]] <- list(inflow = inflow, ps_load = ps_load_list)
      to_do_ids <- setdiff(to_do_ids, id)
    }
  }
  if(length(to_do_ids) == starting_count) {
    warning("Loop stopped: No progress made. Check for network cycles.")
    break
  }
}

print(paste0("Calculated full 'river_map' object size ",
             format(object.size(river_map), units = "MB")))
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 8) Preparing the data for presentation ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for(id in names(river_map)){

  # 1. Update Inflow (Stays the same as it's a single dataframe)
  river_map[[id]]$inflow <- river_map[[id]]$inflow |>
    mutate(area        = round(area/1000000, 1),
           added_area  = round(added_area/1000000, 1),
           inflow_area = round(inflow_area/1000000, 1),
           len         = round(len/1000, 1),
           r_retN      = round(r_retN, 3),
           r_retP      = round(r_retP, 3),
           l_retN      = round(l_retN, 3),
           l_retP      = round(l_retP, 3),
           zero_N_load = round(zero_N_load/1000, 1),
           zero_P_load = round(zero_P_load/1000, 1))
  # 2. Add concentrations and flow information to the inflow dataframe.
  river_map[[id]]$inflow_conc <- df_mod_cl[df_mod_cl$cach_id == id,]

  ps_sums <- sapply(river_map[[id]]$ps_load, function(df_scenario) {
    if(nrow(df_scenario) > 0) {
      return(c(sum_TN_conc = sum(df_scenario$TN_conc_added, na.rm = TRUE),
               sum_TP_conc = sum(df_scenario$TP_conc_added, na.rm = TRUE)))
    } else {
      return(c(sum_TN_conc = 0, sum_TP_conc = 0))
    }
  })

  # Transpose and convert to dataframe to merge/cbind
  ps_sums_df <- as.data.frame(t(ps_sums))

  # Add these as new columns to inflow_conc
  # Note: This assumes inflow_conc has one row per scenario or matches the ps_load names
  river_map[[id]]$inflow_conc <- cbind(river_map[[id]]$inflow_conc, ps_sums_df) |>
    mutate(tn_conc_f = tn_conc + sum_TN_conc,
           tp_conc_f = tp_conc + sum_TP_conc) |>
    mutate(TN_class = case_when(
      tn_conc_f <  2.00  ~ 1,
      tn_conc_f <= 3.00  ~ 2,
      tn_conc_f <= 6.00  ~ 3,
      tn_conc_f <= 12.00 ~ 4,
      tn_conc_f >  12.00 ~ 5,
      TRUE ~ NA_real_),
           TP_class = case_when(
             tp_conc_f <  0.10  ~ 1,
             tp_conc_f <= 0.14  ~ 2,
             tp_conc_f <= 0.23  ~ 3,
             tp_conc_f <= 0.47  ~ 4,
             tp_conc_f >  0.47  ~ 5,
             TRUE ~ NA_real_)) |>
    tibble::rownames_to_column(var = "scenario")

  # 3. Update PS Loads for ALL scenarios
  # We use lapply to loop through scenarios
  river_map[[id]]$ps_load <- lapply(river_map[[id]]$ps_load, function(df_scenario) {

    # Check if the scenario dataframe has data to avoid errors on empty sets
    if(nrow(df_scenario) > 0) {
      ps_n_sum <- sum(df_scenario$TN, na.rm = TRUE)
      ps_p_sum <- sum(df_scenario$TP, na.rm = TRUE)

      # Handle potential division by zero if sum is 0
      df_scenario <- df_scenario |>
        mutate(TN_proc = if(ps_n_sum > 0) round(100*(TN/ps_n_sum), 1) else 0,
               TP_proc = if(ps_p_sum > 0) round(100*(TP/ps_p_sum), 1) else 0,
               TN = round(TN/1000, 3),
               TP = round(TP/1000, 3),
               TN_conc_added = round(TN_conc_added, 3),
               TP_conc_added = round(TP_conc_added, 3))
    }
    return(df_scenario)
  })
}

print(paste0("Final full 'river_map' object size ",
             format(object.size(river_map), units = "MB")))

# ## Test single
# x <- river_map[["8274"]]

## Optionally save the final river_map to an RDS file for use in the Shiny app or other analyses
if(save_river_map){
  saveRDS(river_map, file = file.path(temp_folder, "river_map.rds"))
  message("Final river_map saved to: ", file.path(temp_folder, "river_map.rds"))
}

## Print a summary of the final river_map structure
print("Point source load allocation script finished successfully.")
print(paste0("Final river_map contains ", length(river_map), " catchments."))
