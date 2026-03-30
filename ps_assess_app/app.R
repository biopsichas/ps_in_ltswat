## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 1) Loading required libraries------------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(shiny)
library(leaflet)
library(sf)
library(tidyverse)
library(DT)
library(openxlsx)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 2) Reading data------------------------------------------------
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Data is saved after running ps_assess.R and saved in temp_folder
temp_folder <- "Temp"

## Data path
data_path <- "Data"

## Loading saved objects needed to the app
river_map <- readRDS(file.path("..", temp_folder, "river_map.rds"))

## GIS data for the leaflet maplet
basins <- readRDS(file.path("..", temp_folder, "basins_simp.rds"))
segments <- readRDS(file.path("..", temp_folder, "segments_simp.rds")) |>
  mutate(wb_code = ifelse(!is.na(wbriver_co), paste0("LT", as.character(wbriver_co)), wlake_wb))

## Water body to cach_id relationship, and problematic water bodies that need
## to be checked manually
wb_rel <- file.path("..", data_path, "WB_representative_SWAT_segments_final.xlsx")

## Reading the relationship between water bodies and catchment IDs, as well as a list
## of problematic water bodies that need to be checked manually.
wb_to_cach_id <- read.xlsx(wb_rel, sheet = "cach_id_to_wb")
wb_problematic <- read.xlsx(wb_rel, sheet = "problem_wb")

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 3) Reading data -----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Extract unique water body IDs for the dropdown choices
problem_wb_ids <- unique(wb_problematic$wb_code)

## PS scenarios
ps_scenarios <- names(river_map[[1]]$ps_load)

## Pollutants for dropdown
pollutant <- c("Bendrasis azotas", "Bendrasis fosforas")

ui <- fluidPage(
  titlePanel("Taškinių šaltinių poveikio vertinimas"),

  sidebarLayout(
    # Convert names to numeric, sort them, then convert back to character for the input
    sidebarPanel(
      selectizeInput("select_cach", "Probleminiai vandens telkiniai:",
                     choices  = problem_wb_ids,
                     selected = character(0),
                     options  = list(placeholder = "Pasirinkite VTK...")),
      hr(),
      selectizeInput("select_ps_scenario", "Taškinių šaltinių scenarijai:",
                     choices  = ps_scenarios,
                     selected = ps_scenarios[1],
                     options  = list(placeholder = "Pasirinkite scenarijų...")),
      hr(),
      selectizeInput("select_nutrient", "Teršalas:",
                     choices  = pollutant,
                     selected = pollutant[1],
                     options  = list(placeholder = "Pasirinkite teršalą...")),
      hr()
    ),
    mainPanel(
      leafletOutput("map", height = 500),
      hr(),
      tabsetPanel(
        tabPanel("Baseino informacija", DTOutput("inflow_table")),
        tabPanel("Taškinių krūvių informacija", DTOutput("ps_table"))
      )
    )
  )
)

server <- function(input, output, session) {

  # 1. Create a reactive variable to track which ID the tables should show
  id_for_tables <- reactiveVal("")

  # 2. Update dropdown choices on start
  updateSelectizeInput(session, "select_cach", choices = problem_wb_ids, server = TRUE)

  # 3. When the DROPDOWN changes: look up cach_id from wb_code, update id_for_tables
  observeEvent(input$select_cach, {
    req(input$select_cach != "")

    cach_id <- wb_to_cach_id %>%
      filter(wb_code == input$select_cach) %>%
      pull(cach_id)

    req(length(cach_id) > 0)
    id_for_tables(as.character(cach_id[1]))
  })

  # 4. Reactive for Tables (Listens to id_for_tables — now a cach_id)
  table_data <- reactive({
    req(id_for_tables() %in% names(river_map))
    river_map[[id_for_tables()]]
  })

  # 5. Reactive for PS table filtered by scenario
  ps_data <- reactive({
    req(table_data())
    req(input$select_ps_scenario != "")

    table_data()$ps_load[[input$select_ps_scenario]]
  })

  # 5b. Reactive for PS spatial data
  ps_data_sf <- reactive({
    req(ps_data())

    ps_data() %>%
      select(c(2:8)) %>%
      filter(!is.na(X), !is.na(Y)) %>%
      st_as_sf(coords = c("X", "Y"), crs = 3346) %>%  # LKS94 - adjust if different
      st_transform(4326)
  })

  # 6. Reactive for Map (Listens ONLY to the dropdown)
  # This ensures the blue network only changes when the search box is used
  geo_data <- reactive({
    req(id_for_tables())
    req(id_for_tables() %in% names(river_map))
    relevant_ids <- river_map[[id_for_tables()]]$inflow$cach_id

    # Extract inflow_conc for all relevant catchments and filter by scenario
    inflow_conc <- bind_rows(
      lapply(as.character(relevant_ids), function(id) {
        if (id %in% names(river_map)) river_map[[id]]$inflow_conc
        else NULL
      })
    ) |> filter(scenario == input$select_ps_scenario)

    segments_colored <- segments |>
      filter(id %in% relevant_ids) |>
      left_join(inflow_conc, by = c("id" = "cach_id"))

    list(
      basins   = basins |> filter(cach_id %in% relevant_ids),
      segments = segments_colored,
      ps       = ps_data_sf()
    )
  })

  # 7. Static Base Map
  output$map <- renderLeaflet({
    bbox_init <- st_bbox(segments)
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      fitBounds(
        lng1 = as.numeric(bbox_init[1]), lat1 = as.numeric(bbox_init[2]),
        lng2 = as.numeric(bbox_init[3]), lat2 = as.numeric(bbox_init[4])
      ) %>%
      addPolylines(
        data = segments, group = "Upės", layerId = ~id,
        weight = 1, color = "#34495e", opacity = 0.4,
        label = ~paste("VTK:", wb_code),
        highlightOptions = highlightOptions(weight = 5, color = "red", bringToFront = TRUE),
        labelOptions = labelOptions(sticky = TRUE)
      ) |>
      addLayersControl(

        overlayGroups = c("Upės", "Baseinai", "Parinktos", "Taškiniai šaltiniai"),

        options = layersControlOptions(collapsed = FALSE)

      )
  })

  # 8. Map Observer (Triggered ONLY by dropdown)
  observe({
    req(geo_data())
    g <- geo_data()

    pal <- colorFactor(
      palette = c("blue", "green", "yellow", "orange", "red"),
      levels  = c(1, 2, 3, 4, 5),
      na.color = "grey"
    )

    proxy <- leafletProxy("map") %>%
      clearGroup("Baseinai") %>% clearGroup("Parinktos") %>% clearGroup("Taškiniai šaltiniai")

    proxy %>% addPolygons(data = g$basins, group = "Baseinai", layerId = ~cach_id,
                          weight = 1, color = "grey", fillColor = "blue", fillOpacity = 0.1)
    proxy |>  addPolylines(
      data    = g$segments,
      group   = "Parinktos",
      layerId = ~id,
      weight  = 3,
      color   = ~pal(TP_class),
      opacity = 0.9,
      popup   = ~paste0("<b>Segment:</b> ", id, "<br>",
                        "<b>TN class:</b> ", TN_class, "<br>",
                        "<b>TP class:</b> ", TP_class)
    )

    # Point Sources
    if (nrow(g$ps) > 0) {

      # Scale radius: sqrt scaling, clamped between 4 and 20
      flow <- g$ps[["Nuotėkų kiekis 1000 m3/metus"]]
      radius_scaled <- pmax(4, pmin(20, 4 + sqrt(flow / max(flow, na.rm = TRUE)) * 16))

      proxy %>% addCircleMarkers(
        data = g$ps,
        group = "Taškiniai šaltiniai",
        radius = radius_scaled,
        color       = "black",      # outline
        fillColor   = "red",        # fill
        weight      = 1,            # outline thickness
        stroke = TRUE,
        fillOpacity = 0.8,
        popup = ~paste0(
          "<b>Pavadinimas:</b> ", Pavadinimas, "<br>",
          "<b>Išleistuvo kodas:</b> ", `Išleistuvo kodas`, "<br>",
          "<b>B. azotas (kg/metus):</b> ", `Bendrasis azotas (kg/metus)`, " kg/m<br>",
          "<b>B. fosforas (kg/metus):</b> ", `Bendrasis fosforas (kg/metus)`, " kg/m<br>",
          "<b>Nuotėkų kiekis 1000 m3/metus:</b> ", `Nuotėkų kiekis 1000 m3/metus`, " tūkst. m3"
        )
      )
    }

    # Map fly/zoom only happens here (when dropdown changes)
    bbox <- st_bbox(g$basins)
    proxy %>% flyToBounds(bbox[[1]], bbox[[2]], bbox[[3]], bbox[[4]])
  })

  # 9. Render Tables (Using table_data() which updates on click)
  output$inflow_table <- renderDT({
    req(table_data())
    datatable(table_data()$inflow, options = list(pageLength = 5, scrollX = TRUE))
  })

  # Render PS table using filtered data
  output$ps_table <- renderDT({
    req(ps_data())
    datatable(ps_data(),
              options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })
}

shinyApp(ui, server)
