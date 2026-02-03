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
## 2) Setting parameters to select and data paths-------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Path to the folder where LT SWAT system setups are stored
setup_path <- "E:/LIFE_AAA/swat_lt/Projects/Setup_2020_coarse/Watersheds/"

sc_base <- "basecase"

## Data paths
gis_path <- "Data/GIS/"

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 3) Reading GIS data----------------------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Read basin shapefile
basins <- st_read(dsn = paste0(gis_path,"basin.gdb"), layer = "basin",
                  quiet = T) %>%
  select(c("cach_id","Subbasin","Setup_name", "GRIDCODE"))

## Read reaches shae file
segments <- st_read(paste0(gis_path, "segments_coarse.shp"), quiet = T) %>%
  st_drop_geometry

## Preparing Subbasin Setup_name table
sstable <- basins %>% st_drop_geometry %>% select(Subbasin, Setup_name) %>% unique

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4) Reading SWAT point source input files & calculating loads in parallel-----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cols_to_scale <- c("flo", "sed", "orgn", "sedp", "no3", "solp", "nh3", "no2", "cbod")
tasks <- list()
## Assembling tasts
for (i in seq_len(nrow(sstable))) {
  subb  <- sstable$Subbasin[i]
  setup <- sstable$Setup_name[i]

  base_path <- file.path(setup_path, subb, setup, sc_base)
  recall_rec <- file.path(base_path, "recall.rec")

  if (!file.exists(recall_rec)) next

  dt2 <- SWATreadR:::read_tbl(recall_rec) |>
    filter(str_starts(file, "pt")) |>
    select(name, file)

  dt3 <- SWATreadR:::read_tbl(file.path(base_path, "recall.con")) |>
    filter(str_starts(name, "pt")) |>
    select(name, gis_id, obj_id)

  # Match by name and create task list
  for (j in seq_len(nrow(dt2))) {
    nm <- dt2$name[j]
    tasks <- append(tasks, list(list(
      f_path  = base_path,
      file    = dt2$file[j],
      dt3_row = dt3[dt3$name == nm, ]
    )))
  }
}

# Parallel processing
plan(multisession, workers = 15)
results <- future_lapply(tasks, function(t) {
  read_recall_file(t$f_path, t$file, t$dt3_row)
}, future.seed = TRUE)

# Final assembly
## Nutrients are in kg and sediments in tonnes (in all results)
df_model_in <- bind_rows(results) |>
  left_join(
    basins |> st_drop_geometry() |> select(-GRIDCODE),
    by = c("gis_id" = "cach_id")
  ) |>
  rename(year = yr)


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 5) Reading SWAT+ output files and calculating loads -------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Making task list
tasks <- list()
for (i in seq_len(nrow(sstable))) {
  subb  <- sstable$Subbasin[i]
  setup <- sstable$Setup_name[i]
  if (!file.exists(base_path)) next
  tasks <- append(tasks, list(list(
    subb  = subb,
    setup   = setup
  )))
}

# Parallel processing
plan(multisession, workers = 15)
results <- future_lapply(tasks, function(t) {
  read_recall_output(t$subb, t$setup)
}, future.seed = TRUE)

# Final assembly
df_model_out <- bind_rows(results)

##Cleaning
rm(tasks, results)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 6) Reading apportionment script results -------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Prepare and check the apportionment file path
ap_dir <- file.path(dirname(setup_path), "GIS", "apportionment")
ap_files <- list.files(
  ap_dir,
  pattern = "^apportionment_detailed_basecase_\\d{4}\\.xlsx$",
  full.names = TRUE
)

if (length(ap_files) == 0) {
  stop("No apportionment files found.")
}

ap_meta <- tibble(
  file = ap_files,
  year = str_extract(basename(ap_files), "\\d{4}") |> as.integer()
)


plan(multisession, workers = 10)

df_app <- future_lapply(
  seq_len(nrow(ap_meta)),
  function(i) {
    read_apportionment(
      file = ap_meta$file[i],
      year = ap_meta$year[i]
    )
  },
  future.seed = TRUE
) |> bind_rows()

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 7) Reading LT SWAT point source database tables (post processed)-------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Defining SWAT LT database connection
lt_con <- tryCatch(
  dbConnect(
    RPostgres::Postgres(),
    dbname   = Sys.getenv("LT_DB_NAME"),
    host     = Sys.getenv("LT_DB_HOST"),
    port     = Sys.getenv("LT_DB_PORT"),
    user     = Sys.getenv("LT_DB_USER"),
    password = Sys.getenv("LT_DB_PASSWORD")
  ),
  error = function(e) stop("LT database connection failed: ", e$message)
)
## Check if successful
print(lt_con)

## Load reference tables
ps_catchment <- load_table(lt_con, "point_sources", "ps_catchment")
small_catch  <- load_table(lt_con, "point_sources", "small_catch")

##  Monthly SWAT point source file
ps_monthly_path <- file.path(
  sub(
    "Setup_2020_coarse/Watersheds/?$",
    "Setup_2020_common/Data/PointSources",
    setup_path
  ),
  "monthly.txt"
)

if (!file.exists(ps_monthly_path)) {
  stop("Monthly point source file not found: ", ps_monthly_path)
}

df_raw <- read.table(
  ps_monthly_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)[-1, ]   # drop units row

## Transform SWAT monthly loads
df_ps <- df_raw |>
  mutate(
    DATE = as.Date(DATE, format = "%Y.%m.%d"),
    days_in_month = lubridate::days_in_month(DATE)
  ) |>
  select(-starts_with("rho_")) |>
  mutate(
    across(SS:Discharge, as.numeric),
    across(SS:PORG, ~ (.x * days_in_month) / 1000),  # kg/day → t/month
    SS        = SS / 1000,
    CBOD      = pmax(BOD7 - 4.57 * NH4, 0),
    Discharge = Discharge * days_in_month
  ) |>
  select(-c(BOD7, year, month, days_in_month)) |>
  setNames(c(
    "id", "year", "sed", "no3", "no2", "nh3", "ntot",
    "solp", "ptot", "orgn", "sedp", "flo", "cbod"
  )) |>
  left_join(
    ps_catchment |>
      select(psid, catchmentid) |>
      rename(id = psid, gis_id = catchmentid),
    by = "id"
  )

## Small catchments (DB-based point sources)
df_small_flo <- small_catch[,c("year", "name", "discharge", "y", "x", "catchmentid")] |>
  mutate(coord = paste0(x,"_",y)) |>
  select(-c(x, y)) |>
  group_by(year, coord, name, catchmentid) |>
  summarise(flo = sum(discharge * 1000, na.rm = TRUE), .groups = "drop") |>
  select(-c(coord, name)) |>
  group_by(year, catchmentid) |>
  summarise(flo = sum(flo, na.rm = TRUE), .groups = "drop") |>
  rename(gis_id = catchmentid)  # seconds/year conversion

## Calculating loads from concentrations and discharges for small point sources
df_small <- small_catch |>
  select(year, nutrient, concentration, discharge, catchmentid) |>
  mutate(
    load = concentration * discharge  # seconds/year conversion
  ) |>
  rename(gis_id = catchmentid) |>
  group_by(year, gis_id, nutrient) |>
  summarise(load = sum(load, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(
    id_cols    = c(year, gis_id),
    names_from = nutrient,
    values_from = load,
    values_fill = 0
  ) |>
  mutate(
    CBOD = pmax(BOD7 - 4.57 * NH4, 0),
    SS   = SS / 1000
  ) |>
  select(-BOD7) |>
  setNames(c(
    "year", "gis_id", "nh3", "no2", "no3",
    "ntot", "solp", "ptot", "sed", "cbod"
  )) |>
  mutate(orgn = ntot - (no3 + no2 + nh3),
         sedp = ptot - solp) |>
  left_join(
    df_small_flo,
    by = c("year", "gis_id")
  )


## Combine DB-based point sources
df_raw_in <- df_ps |>
  mutate(year = year(year)) |>
  select(-id) |>
  group_by(year, gis_id) |>
  summarise(
    across(everything(), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  select(names(df_small)) |>
  bind_rows(df_small) |>
  group_by(year, gis_id) |>
  summarise(
    across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

df_ps |>
  mutate(year = year(year)) |>
  select(-id) |>
  group_by(year, gis_id) |>
  summarise(
    across(everything(), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  select(-c(gis_id)) |>
  group_by(year)|>
  summarise(
    across(everything(), \(x) round(sum(x, na.rm = TRUE)/1000, 2)),
    .groups = "drop"
  ) |>
  mutate(id = 1) |>
  select(-year) |>
  group_by(id)|>
  summarise(
    across(everything(), \(x) round(mean(x, na.rm = TRUE), 2)),
    .groups = "drop"
  )

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 8) Processed all data to compare --------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

dt_in <- df_model_in |>
  mutate(ntot = orgn + no3 + no2 + nh3,
         ptot = solp + sedp) |>
  select(year, gis_id, flo, ntot, ptot)

ps_link <-  df_model_in |>
  select(name, gis_id, Subbasin, Setup_name) |>
  mutate(
    name = stringr::str_replace(
      name,
      "^pt(\\d)$",
      "pt0\\1"
    )
  ) |>
  unique()

dt_out <- df_model_out |>
  left_join(
    ps_link,
    by = c("id" = "name", "Subbasin", "Setup_name")
  )|>
  mutate(ntot = orgn + no3 + no2 + nh3,
         ptot = solp + sedp) |>
  select(year, gis_id, flo, ntot, ptot)

dt_ap <- df_app |>
  select(year, gis_id, ntot, ptot)

dt_db <- df_raw_in |>
  mutate(year = as.integer(year),
         gis_id = as.numeric(gis_id)) |>
  select(year, gis_id, flo, ntot, ptot)

## Join all data
dt_all <- dt_in |>
  rename_with(~ paste0(.x, "_in"), -c(year, gis_id)) |>
  full_join(
    dt_out |> rename_with(~ paste0(.x, "_out"), -c(year, gis_id)),
    by = c("year", "gis_id")
  ) |>
  full_join(
    dt_ap |> rename_with(~ paste0(.x, "_ap"), -c(year, gis_id)),
    by = c("year", "gis_id")
  ) |>
  full_join(
    dt_db |> rename_with(~ paste0(.x, "_db"), -c(year, gis_id)),
    by = c("year", "gis_id")
  ) |>
  replace_na(list(
    flo_in = 0, ntot_in = 0, ptot_in = 0,
    flo_out = 0, ntot_out = 0, ptot_out = 0,
    ntot_ap = 0, ptot_ap = 0,
    ntot_db = 0, ptot_db = 0
  )) |>
  select(
    year, gis_id,
    starts_with("flo"),
    starts_with("ntot"),
    starts_with("ptot"))
