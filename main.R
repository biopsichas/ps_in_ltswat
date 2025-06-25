##------------------------------------------------------------------------------
## 1) Loading required libraries
##------------------------------------------------------------------------------

library(tidyverse)
library(sf)
if (!requireNamespace("SWATreadR", quietly = TRUE)) {
  devtools::install_github("chrisschuerz/SWATreadR")
}
##------------------------------------------------------------------------------
## 2) Setting parameters to select and data paths
##------------------------------------------------------------------------------
select_pars <- c("BDS7", "Amonio azotas (NH4-N)", "Bendrasis azotas",
                 "Bendrasis fosforas", "Fosfatinis fosforas (PO4-P)",
                 "Nitratinis azotas (NO3-N)", "Nitritinis azotas (NO2-N)")

## Path to the folder where LT SWAT system setups are stored
setup_path <- "G:/LIFE_AAA/swat_lt/Projects/Setup_2020_coarse/Watersheds/"

## Data paths
data_path <- "Data/"
gis_path <- "Data/GIS/"

## Point source CSV file
point_source_file <- "2023.csv"

## Model setup file name to read
setup_name <- "v_dec_tbl_PAIC9"

## Which channel_sd files to read for yearly - "yr", for annual averages "aa"
channel_sd <- "yr"

##------------------------------------------------------------------------------
## 3) Reading data
##------------------------------------------------------------------------------

# Read point source data
pnt_source_df <- read.csv(paste0(data_path, point_source_file), sep = ";",
                          check.names = FALSE, fileEncoding = "UTF-8")

## Read basin shapefile
basins <- st_read(dsn = paste0(gis_path,"basin.gdb"), layer = "basin",
                  quiet = T) %>%
  select(c("cach_id","Subbasin","Setup_name", "GRIDCODE"))

## Read reaches shae file
segments <- st_read(paste0(gis_path, "segments_coarse.shp"), quiet = T) %>%
  st_drop_geometry

##------------------------------------------------------------------------------
## 4) Cleaning and transforming point source data
##------------------------------------------------------------------------------

# Transform point source data
pnt_df <- pnt_source_df %>%
  select(c("Išleistuvo kodas", "Ūkio subjekto pavadinimas",
           "Išleistuvo koordinatės (LKS'94)", "Nuotekų kiekis, tūkst. m3" ,
           "Laikotarpis", "Teršalo pavadinimas",
           "Vid. metinė (laikotarpio)\u00A0 koncentracija išleidžiamose nuotekose")) %>%
  {
    colnames(.) <- c("code", "name", "coordinates", "wastewater_volume",
                     "period", "pollutant_name", "avg_concentration")
    .
  } %>%
  separate(coordinates, into = c("x", "y"), sep = " ") %>%
  separate(period, into = c("start", "end"), sep = " - ") %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    wastewater_volume = as.numeric(gsub(",", ".", wastewater_volume)),
    avg_concentration = as.numeric(gsub(",", ".", avg_concentration)),
    start = as.Date(start),
    end = as.Date(end)
  ) %>%
  mutate(avg_concentration = case_when(
    pollutant_name == "Amonis (NH4)"     ~ avg_concentration * 14 / 18,
    pollutant_name == "Nitratai (NO3)"   ~ avg_concentration * 14 / 62,
    pollutant_name == "Nitritai (NO2)"   ~ avg_concentration * 14 / 46,
    pollutant_name == "Fosfatai (PO4)"   ~ avg_concentration * 31 / 95,
    TRUE ~ avg_concentration
  )) %>%
  mutate(
    pollutant_name = case_when(
      pollutant_name == "Amonis (NH4)" ~ "Amonio azotas (NH4-N)",
      pollutant_name == "Fosfatai (PO4)" ~ "Fosfatinis fosforas (PO4-P)",
      pollutant_name == "Nitratai (NO3)" ~ "Nitratinis azotas (NO3-N)",
      pollutant_name == "Nitritai (NO2)" ~ "Nitritinis azotas (NO2-N)",
      TRUE ~ pollutant_name  # keep unchanged if not matched
  )) %>%
  filter(pollutant_name %in% select_pars) %>%
  filter(avg_concentration > 0) %>%
  mutate(time_dif = as.numeric(end - start + 1)) %>%
  mutate(load = ifelse(time_dif == max(time_dif, na.rm = TRUE),
                       wastewater_volume * avg_concentration / 1000,
                       (time_dif / max(time_dif, na.rm = TRUE)) *
                         wastewater_volume * avg_concentration / 1000)) %>%
  st_as_sf(coords = c("x", "y"), crs = 3346)

## Joining basins information to point source codes with spatial overlay
pnt_sf <-   st_join(pnt_df["code"] %>% unique, basins, left = TRUE) %>%
  st_drop_geometry %>%
  {
    .[.$code == 1850002, c("cach_id", "Subbasin", "Setup_name", "GRIDCODE")] <- list(21492, "Merkys", "Salcia", 10)
    . <- .[!is.na(.$Subbasin), ]
  }

## Joining point source data with basin information
pst_info <- pnt_df %>%
  st_drop_geometry %>%
  inner_join(pnt_sf, by = "code") %>%
  select(code, name, pollutant_name, load, cach_id, Subbasin, Setup_name, GRIDCODE)

## Preparing Subbasin Setup_name table
sstable <- basins %>% st_drop_geometry %>% select(Subbasin, Setup_name) %>% unique

##------------------------------------------------------------------------------
## 5) Reading SWAT output files and calculating loads
##------------------------------------------------------------------------------

df <- NULL
for(i in 1:dim(sstable)[1]){
  f_path <- paste0(setup_path, sstable[i, "Subbasin"], "/", sstable[i, "Setup_name"],
                   "/", setup_name, "/channel_sd_", channel_sd, ".txt")
  if (file.exists(f_path)) {
    txt <- SWATreadR::read_swat(f_path) %>%
      select(c("unit", "flo_out", "sedp_out", "solp_out", "orgn_out", "nh3_out",
               "no3_out", "no2_out", "cbod_out")) %>%
      group_by(unit) %>%
      summarise(across(everything(), function(x) mean(x, na.rm = TRUE))) %>%
      mutate(ntot_out = orgn_out + nh3_out + no3_out + no2_out,
             ptot_out = sedp_out + solp_out,
             bod7_out = cbod_out * 1.2,
             nh4_out = nh3_out * 1.287) %>%
      pivot_longer(-c(unit, flo_out), names_to = "var_name", values_to = "load_kg_y")%>%
      mutate(conc_mg_l = (load_kg_y * 1e6) / (flo_out * 31536000 * 1000),
             Subbasin =  sstable[i, "Subbasin"],
             Setup_name = sstable[i, "Setup_name"])
    if(!is.null(df)) df <- rbind(df, txt) else df <- txt
  } else {
    if(!grepl("Prieglius", f_path)){
      message(paste0("File ", f_path, " doesn't exist"))
    }
  }
}

##------------------------------------------------------------------------------
## 6) Cleaning up modeling data
##------------------------------------------------------------------------------

## Just to add cach_id info
basins_df <- basins %>%
  st_drop_geometry %>%
  select(cach_id, Subbasin, Setup_name, GRIDCODE)

## Joining SWAT output data with basin information, cleaning names, filtering
df_fix <- df %>%
  mutate(pollutant_name = case_when(
    var_name == "solp_out" ~ "Fosfatinis fosforas (PO4-P)",
    var_name == "no3_out" ~ "Nitratinis azotas (NO3-N)",
    var_name == "nh4_out" ~ "Amonio azotas (NH4-N)",
    var_name == "no2_out" ~ "Nitritinis azotas (NO2-N)",
    var_name == "ntot_out" ~ "Bendrasis azotas",
    var_name == "ptot_out" ~ "Bendrasis fosforas",
    var_name == "bod7_out" ~ "BDS7"
  )) %>%
  filter(!is.na(pollutant_name)) %>%
  rename(GRIDCODE = unit) %>%
  select(pollutant_name, Subbasin, Setup_name, GRIDCODE, load_kg_y) %>%
  left_join(basins_df, by = c("Subbasin", "Setup_name", "GRIDCODE"))

##------------------------------------------------------------------------------
## 7) Get lookup1 table for catch ids that represents the downstream segments of WBs
##------------------------------------------------------------------------------

# Convert to a named list for fast lookup by 'id'
segments_lookup <- segments %>%
  select(id, riverto, wbriver_co, wlake_wb) %>%
  tibble::column_to_rownames("id")

# Prepare an empty results data frame
lookup1 <- tibble(id = numeric(), real_id = numeric())

# Loop over each segment
for (start_id in segments$id) {
  current_id <- start_id
  wb <- wb2 <- "S"  # Initial dummy values to enter the loop

  # Traverse downstream until waterbody changes
  while (TRUE) {
    current_row <- segments_lookup[as.character(current_id), ]

    # Determine current segment's waterbody
    wb <- if (!is.na(current_row$wbriver_co)) current_row$wbriver_co
    else if (!is.na(current_row$wlake_wb)) current_row$wlake_wb
    else wb

    # Get the ID of the downstream segment
    next_id <- current_row$riverto

    # If there is a downstream segment, get its waterbody
    if (next_id != -1 && as.character(next_id) %in% rownames(segments_lookup)) {
      next_row <- segments_lookup[as.character(next_id), ]
      wb2 <- if (!is.na(next_row$wbriver_co)) next_row$wbriver_co
      else if (!is.na(next_row$wlake_wb)) next_row$wlake_wb
      else wb2
    } else {
      wb2 <- "E"  # if end of river system
    }

    # If waterbody changes or flow ends, mark the transition point
    if (wb2 == "E") {
      lookup1 <- add_row(lookup1, id = start_id, real_id = current_id)
      break
    } else if (wb == "S") {
      NULL # Continue, if no waterbody code at the start
    } else if (wb != wb2) {
      lookup1 <- add_row(lookup1, id = start_id, real_id = current_id)
      break
    }
    # Continue traversal
    current_id <- next_id
  }
}

##------------------------------------------------------------------------------
## 8) Calculating final table with percentages of point source impacts
##------------------------------------------------------------------------------

pst_info_prc <- pst_info %>%
  left_join(lookup1, by = c("cach_id" = "id")) %>%
  left_join(df_fix[c("cach_id", "pollutant_name", "load_kg_y")], by = c("pollutant_name", "cach_id")) %>%
  mutate(prc = round(100*load/load_kg_y, 3)) %>%
  left_join(df_fix[c("cach_id", "pollutant_name", "load_kg_y")] %>%
              rename(rload_kg_y = load_kg_y), by = c("pollutant_name", "real_id" = "cach_id")) %>%
  mutate(rprc = ifelse(rload_kg_y > 0, round(100*load/rload_kg_y, 3), prc)) %>%
  left_join(select(segments, id, kadastroid, wbriver_co, wlake_wb), by = c("cach_id" = "id")) %>%
  left_join(select(segments, id, wbriver_co, wlake_wb) %>%
              rename(rwbriver_co = wbriver_co,
                     rwlake_wb = wlake_wb), by = c("real_id" = "id")) %>%
  mutate(wbriver_co = ifelse(is.na(wbriver_co), rwbriver_co, wbriver_co),
         wlake_wb = ifelse(is.na(wlake_wb), rwlake_wb, wlake_wb)) %>%
  select(-c("rwbriver_co", "rwlake_wb"))

colnames(pst_info_prc) <- c("Išleistuvo kodas", "Ūkio subjekto pavadinimas",
                           "Teršalo pavadinimas", "Teršalo kiekis išleidžiamose nuotekose (kg/metus)",
                        "Tiesioginis VT cach_id", "Subbasin", "Setup_name", "GRIDCODE", "Žemiausio taško VT cach_id",
                        "Teršalo kiekis SWAT kg/metus (tiesioginis VT)", "Nuotėkų šaltinių dalis (% tiesioginiame VT)",
                        "Teršalo kiekis SWAT kg/metus (žemiausias VT taškas)", "Nuotėkų šaltinių dalis (%  žemiausias VT taškas)",
                        "Kadastro kodas (tiesioginis VT)", "VT kodas (upė)", "VT kodas (ežeras)")

## Saving the resulting table as csv
write.csv(pst_info_prc, file = paste0("point_source_loads.csv"),
          row.names = FALSE, fileEncoding = "UTF-8", quote = FALSE)

