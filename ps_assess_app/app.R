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

## Load the water body names and types
wb_names <- readRDS(file.path("..", data_path, "wb_names.rds")) |>
  filter(nchar(`VT kodas`) > 6) |>
  mutate(tipas = ifelse(tipas == "U", "Upė", "Ežeras/Tvenkinys"))

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
## 3) Setting page -----
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
      hr(),
      uiOutput("selection_info")
    ),
    mainPanel(
      leafletOutput("map", height = 500),
      hr(),
      tabsetPanel(
        tabPanel("Taškinių krūvių informacija", DTOutput("ps_table")),
        tabPanel("Baseino informacija", DTOutput("inflow_table"))
      )
    )
  )
)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4) Server behavoir -----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
server <- function(input, output, session) {

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4.1) Reactive variable to track selected cach_id for tables ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  #  Create a reactive variable to track which ID the tables should show
  id_for_tables <- reactiveVal("")

  # Update dropdown choices on start
  updateSelectizeInput(session, "select_cach", choices = problem_wb_ids, server = TRUE)

  # When the DROPDOWN changes: look up cach_id from wb_code, update id_for_tables
  observeEvent(input$select_cach, {
    req(input$select_cach != "")

    cach_id <- wb_to_cach_id |>
      filter(wb_code == input$select_cach) |>
      pull(cach_id)

    req(length(cach_id) > 0)
    id_for_tables(as.character(cach_id[1]))
  })

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4.2) Tables and map data ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Reactive for Tables (Listens to id_for_tables — now a cach_id)
  table_data <- reactive({
    req(id_for_tables() %in% names(river_map))
    river_map[[id_for_tables()]]
  })

  # Reactive for Tables (Listens to id_for_tables — now a cach_id)
  table_data_conc <- reactive({
    req(table_data())
    req(input$select_ps_scenario != "")
    df <- table_data()$inflow_conc
    ## drop = FALSE to ensure it stays a data frame even if 1 row
    df[df$scenario == input$select_ps_scenario, , drop = FALSE]
  })

  # Reactive for PS table filtered by scenario
  ps_data <- reactive({
    req(table_data())
    req(input$select_ps_scenario != "")

    table_data()$ps_load[[input$select_ps_scenario]]
  })

  # Reactive for PS spatial data
  ps_data_sf <- reactive({
    d <- ps_data()
    if (is.null(d) || nrow(d) == 0) return(NULL)

    d |>
      select(c(2:8)) |>
      filter(!is.na(X), !is.na(Y)) |>
      st_as_sf(coords = c("X", "Y"), crs = 3346) |>  # LKS94 - adjust if different
      st_transform(4326)
  })

  # Reactive for Map (Listens ONLY to the dropdown)
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


    common_cols <- c("id", "flo_out")

    nutrient_cols <- switch(input$select_nutrient,
                            "Bendrasis azotas"   = c("tn_conc", "sum_TN_conc", "tn_conc_f", "TN_class"),
                            "Bendrasis fosforas" = c("tp_conc", "sum_TP_conc", "tp_conc_f", "TP_class")
    )

    segments_colored <- segments |>
      filter(id %in% relevant_ids) |>
      left_join(inflow_conc, by = c("id" = "cach_id")) |>
      select(all_of(c(common_cols, nutrient_cols)))

    list(
      basins   = basins |> filter(cach_id %in% relevant_ids),
      segments = segments_colored,
      ps       = ps_data_sf()
    )
  })

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4.3) Map ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Static Base Map
  output$map <- renderLeaflet({
    bbox_init <- st_bbox(segments)
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron,
                       group = "CartoDB") |>
      addProviderTiles(providers$Esri.WorldImagery,
                       group = "Ortofoto") |>
      fitBounds(
        lng1 = as.numeric(bbox_init[1]), lat1 = as.numeric(bbox_init[2]),
        lng2 = as.numeric(bbox_init[3]), lat2 = as.numeric(bbox_init[4])
      ) |>
      addPolylines(
        data = segments, group = "Upės", layerId = ~id,
        weight = 1, color = "blue", opacity = 0.4,
        label = ~paste("VTK:", wb_code),
        highlightOptions = highlightOptions(weight = 5, color = "red", bringToFront = TRUE),
        labelOptions = labelOptions(sticky = TRUE)
      ) |>
      addLayersControl(
        baseGroups    = c("CartoDB", "Ortofoto"),
        overlayGroups = c("Upės", "Parinkti baseinai", "Parinktos upės", "Taškiniai šaltiniai"),
        options       = layersControlOptions(collapsed = FALSE)
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
    nutrient_col <- switch(input$select_nutrient,
                           "Bendrasis azotas"   = "TN_class",
                           "Bendrasis fosforas" = "TP_class"
    )


    proxy <- leafletProxy("map") |>
      clearGroup("Parinkti baseinai") |> clearGroup("Parinktos upės") |> clearGroup("Taškiniai šaltiniai")

    proxy |> addPolygons(data = g$basins, group = "Parinkti baseinai", layerId = ~cach_id,
                          weight = 1, color = "grey", fillColor = "blue", fillOpacity = 0.1)

    class_label <- function(class_val) {
      dplyr::case_when(
        class_val == 1 ~ "Labai gera",
        class_val == 2 ~ "Gera",
        class_val == 3 ~ "Vidutinė",
        class_val == 4 ~ "Bloga",
        class_val == 5 ~ "Labai bloga",
        TRUE           ~ "Nežinoma"
      )
    }

    segment_popup <- switch(input$select_nutrient,
                            "Bendrasis azotas" = ~paste0(
                              "<b>Baseinėlio ID:</b> ", id, "<br>",
                              "<b>Debitas:</b> ", round(flo_out, 3), " m3/s<br>",
                              "<b>Bazinė TN konc.:</b> ", round(tn_conc, 2), " N mg/l<br>",
                              "<b>Taškinių pridedama TN konc.:</b> ", round(sum_TN_conc, 2), " N mg/l<br>",
                              "<b>Galutinė TN konc.:</b> ", round(tn_conc_f, 2), " N mg/l<br>",
                              "<b>VT būklė pagal TN:</b> ", class_label(TN_class)
                            ),
                            "Bendrasis fosforas" = ~paste0(
                              "<b>Baseinėlio ID:</b> ", id, "<br>",
                              "<b>Debitas:</b> ", round(flo_out, 3), " m3/s<br>",
                              "<b>Bazinė TP konc.:</b> ", round(tp_conc, 3), " P mg/l<br>",
                              "<b>Taškinių pridedama TP konc.:</b> ", round(sum_TP_conc, 3), " P mg/l<br>",
                              "<b>Galutinė TP konc.:</b> ", round(tp_conc_f, 3), " P mg/l<br>",
                              "<b>VT būklė pagal TP:</b> ", class_label(TP_class)
                            )
    )

    ## To make it easier to click on thin lines, we add a transparent thicker line underneath for popups and highlighting
    proxy |> addPolylines(
      data    = g$segments,
      group   = "Parinktos upės",
      layerId = paste0("buffer_", g$segments$id),
      weight  = 15,
      color   = "black",
      opacity = 0,
      popup   = segment_popup
    )

    proxy |>  addPolylines(
      data    = g$segments,
      group   = "Parinktos upės",
      layerId = ~id,
      weight  = 3,
      color   = ~pal(get(nutrient_col)),
      opacity = 0.9,
      popup   = segment_popup
    )

    proxy |> addLegend(
      pal = pal,
      values = c(1, 2, 3, 4, 5),
      title = "Būklės klasė",
      position = "bottomright",
      layerId = "legend",
      labFormat = labelFormat(
        transform = function(x) {
          # This looks up the text for each numeric value in the legend
          sapply(x, class_label)
        }
      )
    )

    # Point Sources
    if (!is.null(g$ps) && nrow(g$ps) > 0) {

      # Scale radius: sqrt scaling, clamped between 4 and 20
      flow <- g$ps[["Nuotėkų kiekis 1000 m3/metus"]]
      radius_scaled <- pmax(4, pmin(20, 4 + sqrt(flow / max(flow, na.rm = TRUE)) * 16))

      proxy |> addCircleMarkers(
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
          "<b>B. azotas (kg/metus):</b> ", `Bendrasis azotas (kg/metus)`, "<br>",
          "<b>B. fosforas (kg/metus):</b> ", `Bendrasis fosforas (kg/metus)`, "<br>",
          "<b>Nuotėkų kiekis 1000 m3/metus:</b> ", `Nuotėkų kiekis 1000 m3/metus`
        )
      )
    }

    # Map fly/zoom only happens here (when dropdown changes)
    bbox <- st_bbox(g$basins)
    proxy |> flyToBounds(bbox[[1]], bbox[[2]], bbox[[3]], bbox[[4]])

  })

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4.4) Info table  ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  output$selection_info <- renderUI({
    req(input$select_cach != "")

    wb_info <- wb_names %>%
      filter(`VT kodas` == input$select_cach)

    wb_name <- if (nrow(wb_info) > 0) wb_info$Pavadinimas[1] else "—"
    wb_type <- if (nrow(wb_info) > 0) wb_info$tipas[1]       else "—"

    prob_info <- wb_problematic %>%
      filter(wb_code == input$select_cach)

    imp_tn <- if (nrow(prob_info) > 0) prob_info$imp_tn[1] else NA
    imp_tp <- if (nrow(prob_info) > 0) prob_info$imp_tp[1] else NA

    tn_badge <- if (isTRUE(imp_tn))
      tags$span("⚠ TN problema", style = "background:#fde8e8; color:#c0392b; padding:1px 7px; border-radius:10px; font-size:14px; font-weight:600;")
    else
      tags$span("✓ TN gerai",    style = "background:#eafaf1; color:#27ae60; padding:1px 7px; border-radius:10px; font-size:14px; font-weight:600;")

    tp_badge <- if (isTRUE(imp_tp))
      tags$span("⚠ TP problema", style = "background:#fde8e8; color:#c0392b; padding:1px 7px; border-radius:10px; font-size:14px; font-weight:600;")
    else
      tags$span("✓ TP gerai",    style = "background:#eafaf1; color:#27ae60; padding:1px 7px; border-radius:10px; font-size:14px; font-weight:600;")

    totals <- table_data()$inflow |>
      summarise(
        total_area = sum(area, na.rm = TRUE),
        added_area = sum(added_area, na.rm = TRUE),
        inflow_area = sum(inflow_area, na.rm = TRUE),
        len = sum(len, na.rm = TRUE)
      )

    totals_id <- table_data()$inflow[table_data()$inflow$cach_id == id_for_tables(),]

    # Scenario-Specific Concentrations
    conc_data <- table_data_conc()
    current_flow <- conc_data$flo_out[1]
    seconds_in_year <- 365.25 * 24 * 3600
    ## Calculate Maximum Possible Load (t/metus) if concentrations were at the threshold level for the entire flow
    max_load_n <- (current_flow * 3 * seconds_in_year) / 1e6
    max_load_p <- (current_flow * 0.14 * seconds_in_year) / 1e6

    # Calculate Current Yearly Load (t/metus) for comparison
    current_load_n <- (current_flow * conc_data$tn_conc_f[1] * seconds_in_year) / 1e6
    current_load_p <- (current_flow * conc_data$tp_conc_f[1] * seconds_in_year) / 1e6

    # Required reduction to meet thresholds
    req_red_n <- max_load_n - current_load_n
    req_red_p <- max_load_p - current_load_p

    # Required reduction to meet thresholds set to 0 if above 0
    req_red_n <- if(req_red_n < 0) req_red_n else 0
    req_red_p <- if(req_red_p < 0) req_red_p else 0

    # Calculate the contribution of point sources to the current load
    ps_load_n <- (current_flow * conc_data$sum_TN_conc[1] * seconds_in_year) / 1e6
    ps_load_p <- (current_flow * conc_data$sum_TP_conc[1] * seconds_in_year) / 1e6

    # Calculate the required percentage reduction in point source loads to meet the thresholds, only if reduction is needed and possible
    # For TN
    if(req_red_n == 0) {
      tn_requirement <- span(style = "color:#276940; font-weight:bold; font-size:14px;", "mažinimas nereikalingas!")
    } else if (abs(req_red_n) <= ps_load_n) {
      tn_requirement <- round(100 * req_red_n/ ps_load_n, 1)
      tn_requirement <- span(style = "color:#D4660F; font-weight:bold; font-size:14px;", paste(tn_requirement, "%"))
    } else {
      tn_requirement <- span(style = "color:#C41800; font-weight:bold; font-size:14px;",  "neįmanomas!!!")
    }

    # For TP
    if(req_red_p == 0) {
      tp_requirement <- span(style = "color:#276940; font-weight:bold; font-size:14px;", "mažinimas nereikalingas!")
    } else if (abs(req_red_p) <= ps_load_p) {
      tp_requirement <- round(100 * req_red_p/ ps_load_p, 1)
      tp_requirement <- span(style = "color:#D4660F; font-weight:bold; font-size:14px;", paste(tp_requirement, "%"))
    } else {
      tp_requirement <- span(style = "color:#C41800; font-weight:bold; font-size:14px;",  "neįmanomas!!!")
    }

    ps_number <- if (!is.null(ps_data()) && is.data.frame(ps_data())) {
      nrow(ps_data())
    } else {
      0
    }

    div(
      style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;
             padding:8px 14px; font-size:13px; line-height:1.9;",
      div(style = "font-size:14px; margin-bottom:8px; color:#2c3e50;",
          tags$b("VANDENS TELKINIO INFORMACIJA")),
      tags$hr(style="margin: 5px 0 10px 0; border-top: 1px solid #ccc;"),
      tags$b("Pavadinimas: "), wb_name, tags$br(),
      tags$b("Tipas: "), wb_type, tags$br(),
      tags$b("VT kodas: "), input$select_cach, tags$br(),
      tags$b("Kas negerai?: "), tn_badge, " ", tp_badge,  tags$br(),
      tags$b("Plotas LT: "), totals$total_area,  " km2", tags$br(),
      tags$b("Visas plotas: "), totals$total_area + totals$added_area + totals$inflow_area, " km2", tags$br(),
      tags$b("Segmentų ilgis: "), totals$len, " km", tags$br(),
      tags$b("Taškinių šaltinių sk.: "), ps_number, tags$br(),
      tags$b("Bazinis sc. NT krūvis: "), totals_id$zero_N_load, " t/m", tags$br(),
      tags$b("Bazinis sc. PT krūvis:"), totals_id$zero_P_load, " t/m", tags$br(),
      tags$b("Debitas: "), round(current_flow, 3), " m3/s", tags$br(),
      tags$hr(style="margin: 5px 0; border-top: 1px dashed #ccc;"),
      # Nitrogen Info
      div(style = "margin-top:5px; color:#2c3e50; font-weight:bold;", "B. azotas:"),
      tags$ul(style = "margin-bottom:5px; padding-left:20px;",
              tags$li("Esama konc.: ", round(conc_data$tn_conc_f[1], 2), " N mg/l"),
              tags$li("Maks. leistinas krūvis: ", tags$b(round(max_load_n, 2)), " t/m"),
              tags$li("Esamas metinis krūvis: ", round(current_load_n, 2), " t/m"),
              tags$li("Reikalingas sumažinimas: ", round(req_red_n, 2), " t/m"),
              tags$li("Taškinių šaltinių indėlis: ", round(ps_load_n, 2), " t/m"),
              tags$li("Reikalingas min. šaltinių sumažinimas pasiekti g.b.:",
                      tags$br(),
                      tn_requirement
              )
      ),

      # Phosphorus Info
      div(style = "color:#2c3e50; font-weight:bold;", "B. fosforas:"),
      tags$ul(style = "margin-bottom:5px; padding-left:20px;",
              tags$li("Esama konc.: ", round(conc_data$tp_conc_f[1], 3), " P mg/l"),
              tags$li("Maks. leistinas krūvis: ", tags$b(round(max_load_p, 3)), " t/m"),
              tags$li("Esamas metinis krūvis: ", round(current_load_p, 3), " t/m"),
              tags$li("Reikalingas sumažinimas: ", round(req_red_p, 2), " t/m"),
              tags$li("Taškinių šaltinių indėlis ", round(ps_load_p, 3), " t/m"),
              tags$li("Reikalingas min. šaltinių sumažinimas pasiekti g.b.:",
                      tags$br(),
                      tp_requirement
              )
      ),
      tags$hr(style="margin: 5px 0; border-top: 1px dashed #ccc;"),

    )
  })

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4.5) Point source table  ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Render Tables (Using table_data() which updates on click)
  # Render PS table using filtered data
  output$ps_table <- renderDT({
    req(ps_data())

    if (nrow(ps_data()) == 0) {
      return(datatable(data.frame(Pranešimas = "Šiame scenarijuje taškinių šaltinių nėra"),
                       rownames = FALSE, options = list(dom = 't')))
    }

    common_cols <- c("Išleistuvo kodas", "Pavadinimas",
                     "Nuotėkų kiekis 1000 m3/metus")

    nutrient_cols <- switch(input$select_nutrient,
                            "Bendrasis azotas"   = c("Bendrasis azotas (kg/metus)", "TN", "TN_conc_added"),
                            "Bendrasis fosforas" = c("Bendrasis fosforas (kg/metus)", "TP", "TP_conc_added")
    )

    rename_map <- c(
      "Galutiniame taške (N kg/metus)" = "TN",
      "Pridedama konc. (N mg/l)"       = "TN_conc_added",
      "Galutiniame taške (P kg/metus)" = "TP",
      "Pridedama konc. (P mg/l)"       = "TP_conc_added"
    )

    ps_data() |>
      mutate(
        TN = TN * 1000,
        TP = TP * 1000
      ) |>
      select(all_of(c(common_cols, nutrient_cols))) |>
      rename(any_of(rename_map)) |>
      arrange(desc(across(last_col()))) |>
      datatable(
        extensions = 'Buttons', # 1. Enable Buttons extension
        rownames = FALSE,
        options = list(
          pageLength = 8,
          scrollX = TRUE,
          dom = 'Bfrtip',       # 2. Add 'B' to the layout to show the buttons
          buttons = list(       # 3. Define which buttons to show
            list(
              extend = 'excel',
              filename = paste0("Taskiniai_saltiniai_", input$select_cach)
            ),
            list(
              extend = 'csv',
              filename = paste0("Taskiniai_saltiniai_", input$select_cach)
            ),
            'pdf'
          )
        )
      )
  })

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 4.6) Basin table  ----
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  output$inflow_table <- renderDT({
    req(table_data())

    common_cols <- c("cach_id", "fraction", "area", "added_area", "inflow_area", "len")

    nutrient_cols <- switch(input$select_nutrient,
                            "Bendrasis azotas"   = c("r_retN", "l_retN", "zero_N_load"),
                            "Bendrasis fosforas" = c("r_retP", "l_retP", "zero_P_load")
    )

    rename_map <- c(
      "Baseinėlio ID" = "cach_id",
      "Dalis skirta VT" = "fraction",
      "Baseinėlio plotas (km2)" = "area",
      "Pridėtas plotas (km2)" = "added_area",
      "Įtekantis plotas (km2)" = "inflow_area",
      "Atkarpos ilgis (km)" = "len",
      "N ret. koef. upėje" = "r_retN",
      "N ret. koef. ežere" = "l_retN",
      "P ret. koef. upėje" = "r_retP",
      "P ret. koef. ežere" = "l_retP",
      "N krūvis be PS (kg/metus)" = "zero_N_load",
      "P krūvis be PS (kg/metus)" = "zero_P_load"
    )
    table_data()$inflow |>
      select(any_of(c(common_cols, nutrient_cols))) |>
      rename(any_of(rename_map)) |>
      datatable(
        extensions = 'Buttons', # Enable export buttons
        options = list(
          pageLength = 5,
          scrollX = TRUE,
          dom = 'Bfrtip',       # Layout: Buttons, filter, processing, table, info, pagination
          buttons = list(
            list(
              extend = 'excel',
              filename = paste0("Baseino_duomenys_", input$select_cach)
            ),
            list(
              extend = 'csv',
              filename = paste0("Baseino_duomenys_", input$select_cach)
            ),
            'pdf'
          )
        )
      )
  })
}

shinyApp(ui, server)
