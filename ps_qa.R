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
setup_path <- "F:/Setup_20220708/Projects/Setup_2020_coarse/Watersheds/"

sc_base <- "WFD_pressures_2020_2024_base"

## Data paths
gis_path <- "Data/GIS/"

## Apportionment file name base, without the year part (_YYYY.xlsx)
ap_name <- "apportionment_detailed_WFD_pressures_2020_2024_base"

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
plan(multisession, workers = 10)
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

### Preparing maximum point source could values table needed for reading outputs
## to avoid reading extint file outputs
ps_max <- df_model_in |>
  select(name, Subbasin, Setup_name) |>
  mutate(val = readr::parse_number(name)) |>
  group_by(Subbasin, Setup_name) |>
  summarise(
    # Specifically target 'val' or columns that are actually numeric
    across(where(is.numeric), \(x) max(x, na.rm = TRUE)),
    .groups = "drop"
  )

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
plan(multisession, workers = 10)
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
  pattern = paste0("^", ap_name, "_", "\\d{4}\\.xlsx$"),
  full.names = TRUE
)

if (length(ap_files) == 0) {
  stop("No apportionment files found.")
}

ap_meta <- tibble(
  file = ap_files,
  year = str_extract(basename(ap_files), "\\d{4}(?=\\D*$)") |> as.integer()
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
    user     = Sys.getenv("LT_DB_USER")
  ),
  error = function(e) stop("LT database connection failed: ", e$message)
)
## Check if successful
print(lt_con)

## Load reference tables
ps_catchment <- load_table(lt_con, "point_sources", "ps_catchment")
small_catch  <- load_table(lt_con, "point_sources", "small_up2024")

## Values taken from pointsource_utils.py
fill_gaps <- data.frame(type = c(700,600,500,300,313,312,311,307,306,305,304,303,
                                 302,101,100,0,301,200,900,400,3130,3131),
                        Ntot = c(1.2,1.2,0.8,1.9,1.7,2.0,2.1,1.9,1.2,2.2,2.2,
                                   1.5,1.6,0.97,0.97,0.8,1.9,1.9,1.2,0.8,1.7,1.9),
                        NH4 = c(0.4,0.1,0.8,0.7,0.25,0.4,0.5,0.7,0.6,0.5,0.4,
                                  0.6,0.6,0.6,0.6,0.4,0.7,0.7,0.4,0.8,0.25,0.7),
                        NO3 = c(0.5,0.8,0.02,0.2,0.4,0.4,0.3,0.2,0.1,0.3,0.3,
                                  0.13,0.15,0.1,0.1,0.3,0.2,0.2,0.5,0.02,0.4,0.2),
                        Ptot = c(0.3,0.3,0.2,0.16,0.2,0.5,0.4,0.4,0.3,0.4,0.4,
                                   0.2,0.4,0.14,0.14,0.1,0.16,0.15,0.3,0.2,0.16,0.2))|>
  mutate(PO4 = Ptot * 0.7) |>
  pivot_longer(cols = -type, names_to = "nutrient",
               values_to = "concentration_base")

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
    across(SS:PORG, ~ (.x * days_in_month) / 1000),  # g/day → kg/month
    SS        = SS / 1000, # kg/month to t/month
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
  )|>
  mutate(across(where(bit64::is.integer64), as.numeric))

## Add cach_id
small_catch <- small_catch |>
  select(year, name, type, nutrient, concentration, discharge, x, y) %>%
  st_as_sf(coords = c("x", "y"), crs = 3346, remove = FALSE) |>
  st_join(basins[,c("cach_id")], left = TRUE) |>
  mutate(type = as.integer(type)) |>
  st_drop_geometry()

## Small catchments (DB-based point sources)
df_small_flo <- small_catch[,c("year", "name", "type", "discharge", "y", "x", "cach_id")] |>
  mutate(coord = paste0(x,"_",y)) |>
  select(-c(x, y)) |>
  group_by(year, coord, name, cach_id, type, discharge) |>
  unique()|>
  summarise(flo = sum(discharge * 1000, na.rm = TRUE), .groups = "drop") |> ##
  select(-c(coord, name)) |>
  group_by(year, cach_id) |>
  summarise(flo = sum(flo, na.rm = TRUE), .groups = "drop") |>
  rename(gis_id = cach_id)


## Calculating loads from concentrations and discharges for small point sources
df_small <- small_catch |>
  mutate(coord = paste0(x,"_",y)) |>
  select(name, coord, year, type, nutrient, concentration, discharge, cach_id) |>
  left_join(fill_gaps, by = c("type", "nutrient")) |>
  # Fill gaps only if concentration is NA or zero and there is a base value to fill with
  mutate(concentration = ifelse((is.na(concentration) | concentration <= 0) &
                                  !is.na(concentration_base) , concentration_base, concentration)) |>
  select(-concentration_base) |>
  mutate(concentration = ifelse(is.na(concentration) | concentration < 0, 0, concentration))

## First aggregate by name, year, type, coord, cach_id, nutrient to get weighted mean concentration and total discharge
df_agg <- df_small |>
  group_by(name, year, type, coord, cach_id, nutrient) |>
  summarise(
    concentration = weighted.mean(concentration, discharge, na.rm = TRUE),
    discharge = sum(discharge, na.rm = TRUE),
    type = first(type),
    .groups = "drop"
  ) |>
  # Reshape to wide format to apply nutrient-specific rules
  pivot_wider(
    id_cols    = c(name, coord, year, type, cach_id, discharge),
    names_from = nutrient,
    values_from = concentration,
    values_fill = 0
  ) |>
  mutate(PO4 = ifelse(PO4>Ptot, 0, PO4),
         NO3 = ifelse(NO3>Ntot, 0, NO3),
         NH4 = ifelse(NH4>Ntot, 0, NH4),
         NO2 = ifelse(NO2>Ntot, 0, NO2),
         NO3 = ifelse(NO3>Ntot - NH4 - NO2, 0, NO3),
         NH4 = ifelse(NH4>Ntot - NO3 - NO2, 0, NH4),
         NO2 = ifelse(NO2>Ntot - NO3 - NH4, 0, NO2)) |>
  pivot_longer(cols = -c(name, coord, year, type, cach_id, discharge),
               names_to = "nutrient",
               values_to = "concentration") |>
  ## Filling gaps again for all variables, which have not been in the database.
  left_join(fill_gaps, by = c("type", "nutrient")) |>
  mutate(concentration = ifelse(concentration == 0, concentration_base, concentration)) |>
  select(-concentration_base) |>
  mutate(concentration = ifelse(is.na(concentration) | concentration < 0, 0, concentration))

## Preparing loads for small point sources, applying nutrient-specific rules and filling gaps as needed
df_small_filled <- df_agg |>
  select(year, nutrient, concentration, discharge, cach_id) |>
  mutate(
    load = concentration * discharge  # seconds/year conversion
  ) |>
  rename(gis_id = cach_id) |>
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
  mutate(gis_id = as.numeric(gis_id)) |>
  bind_rows(df_small) |>
  group_by(year, gis_id) |>
  summarise(
    across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

# ## Balance check function
print("Big point sources")
check_balance(df_ps)
print("Small point sources")
check_balance(df_small)
print("Joinned data")
check_balance(df_raw_in)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 8) Processed all data to compare --------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

dt_in <- df_model_in |>
  mutate(ntot = orgn + no3 + no2 + nh3,
         ptot = solp + sedp) |>
  select(year, gis_id, flo, ntot, ptot)

## Create point source name to gis_id link table
ps_link <-  df_model_in |>
  select(name, gis_id, Subbasin, Setup_name) |>
  mutate(id = readr::parse_number(name)) |>
  select(-name) |>
  unique()

dt_out <- df_model_out |>
  left_join(ps_link,
            by = c("id", "Subbasin", "Setup_name")
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
    flo_db = 0, ntot_db = 0, ptot_db = 0
  )) |>
  select(
    year, gis_id,
    starts_with("flo"),
    starts_with("ntot"),
    starts_with("ptot"))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 9) Quality assurance checks -------------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Print the final sums
yr_sel <- seq(1997, 2019)
print("Model input (*.rec files)")
check_balance(filter(dt_in, year %in% yr_sel))
print("Model output (recall_yr.txt files)")
check_balance(filter(dt_out, year %in% yr_sel))
print("Apportionment results (apportionment_detailed_basecase_YYYY.xlsx files)")
check_balance(filter(dt_ap, year %in% yr_sel))
print("Raw input data (monthly.txt and PostGress database)")
check_balance(filter(dt_db, year %in% yr_sel))

## Calculate relative differences and filter for significant discrepancies
x <- dt_all |> select(ends_with(c("gis_id", "year", "_in", "_db"))) |>
  select(
    year, gis_id,
    starts_with("flo"),
    starts_with("ntot"),
    starts_with("ptot")) |>
  mutate(q = round(1 - (flo_in/flo_db), 2),
         q = ifelse(is.infinite(q) | is.na(q), 0, q),
         n = round(1 - (ntot_in/ntot_db), 2),
         n = ifelse(is.infinite(n) | is.na(n), 0, n),
         p = round(1 - (ptot_in/ptot_db), 2),
         p = ifelse(is.infinite(p) | is.na(p), 0, p))


xx <- x |>
  filter(year < 2020) |>
  filter(q != 0 | n != 0 | p != 0) |>
  filter(q < 0.9 | n < 0.9 | p < 0.9)
