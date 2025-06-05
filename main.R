library(tidyverse)
library(sf)
select_pars <- c("BDS7", "Amonio azotas (NH4-N)", "Bendrasis azotas", "Bendrasis fosforas", "Fosfatinis fosforas (PO4-P)", "Nitratinis azotas (NO3-N)", "Nitritinis azotas (NO2-N)")

x <- read.csv("2023.csv", sep = ";", check.names = FALSE, fileEncoding = "UTF-8") %>%
  select(c("Išleistuvo kodas", "Ukio subjekto pavadinimas", "Išleistuvo koordinatės (LKS'94)", "Nuotekų kiekis, tūkst. m3" , "Laikotarpis", "Teršalo pavadinimas", "Vid. metinė (laikotarpio)\u00A0 koncentracija išleidžiamose nuotekose")) %>%
  {
    colnames(.) <- c("code", "name", "coordinates", "wastewater_volume", "period", "pollutant_name", "avg_concentration")
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
                       (time_dif / max(time_dif, na.rm = TRUE)) * wastewater_volume * avg_concentration / 1000)) %>%
  st_as_sf(coords = c("x", "y"), crs = 3346) %>%
  st_transform(4326)


y <- st_read("basin.shp", quiet = T) %>%
  st_transform(4326) %>%
  select(c("cach_id","Subbasin","Setup_name", "GRIDCODE"))

points <- st_join(x["code"] %>% unique, y, left = TRUE) %>%
  st_drop_geometry

points[points$code == 1850002, c("cach_id", "Subbasin", "Setup_name", "GRIDCODE")] <- list(21492, "Merkys", "Salcia", 10)
points <- points[!is.na(points$Subbasin),]

pst_info <- x %>%
  st_drop_geometry %>%
  inner_join(points, by = "code") %>%
  select(code, name, pollutant_name, load, cach_id, Subbasin, Setup_name, GRIDCODE)


v <- y %>% st_drop_geometry %>% select(Subbasin, Setup_name) %>% unique

setup_path <- "G:/LIFE_AAA/swat_lt/Projects/Setup_2020_coarse/Watersheds/"

df <- NULL
for(i in 1:dim(v)[1]){
  f_path <- paste0(setup_path, v[i, "Subbasin"], "/", v[i, "Setup_name"], "/v_dec_tbl_PAIC9/channel_sd_aa.txt")
  if (file.exists(f_path)) {
    txt <- SWATreadR::read_swat(f_path) %>%
      select(c("unit", "flo_out", "sedp_out", "solp_out", "orgn_out", "nh3_out", "no3_out", "no2_out", "cbod_out"))%>%
      mutate(ntot_out = orgn_out + nh3_out + no3_out + no2_out,
             ptot_out = sedp_out + solp_out,
             bod7_out = cbod_out * 1.2,
             nh4_out = nh3_out * 1.287) %>%
      pivot_longer(-c(unit, flo_out), names_to = "var_name", values_to = "load_kg_y")%>%
      mutate(conc_mg_l = (load_kg_y * 1e6) / (flo_out * 31536000 * 1000),
             Subbasin =  v[i, "Subbasin"],
             Setup_name = v[i, "Setup_name"])
    if(!is.null(df)) df <- rbind(df, txt) else df <- txt
  } else {
    message(paste0("File ", f_path, " doesn't exist"))
  }
}


df0 <- df %>%
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
  select(pollutant_name, Subbasin, Setup_name, GRIDCODE, load_kg_y)


segments <- st_read("segments_coarse.shp", quiet = T) %>%
  st_drop_geometry %>%
  select(id, kadastroid, wbriver_co, wlake_wb)

pst_info_prc <- pst_info %>%
  left_join(df0, by = c("pollutant_name", "Subbasin", "Setup_name", "GRIDCODE")) %>%
  mutate(prc = round(100*load/load_kg_y, 3)) %>%
  left_join(segments, by = c("cach_id" = "id"))





ggplot(txt, aes(x = 0, y = conc_mg_l)) +
  geom_boxplot(outlier.shape = NA) +  # hide outliers (optional)
  facet_wrap(~ var_name, scales = "free_y") +  # independent y-axis per facet
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),       # remove x-axis text
    axis.ticks.x = element_blank(),      # remove x-axis ticks
    axis.title.x = element_blank()       # remove x-axis title
  )



