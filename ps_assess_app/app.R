library(shiny)
library(leaflet)
library(sf)
library(tidyverse)
library(DT)

# --- Assumption: You have these objects in your environment ---

# 1. Get the names
all_ids <- names(river_map)

# 2. Sort them numerically (handles 1, 10, 100, 1000000 correctly)
sorted_ids <- stringr::str_sort(all_ids, numeric = TRUE)

# 3. Create the choice list
choices_list <- c("Select an ID" = "", sorted_ids)

ui <- fluidPage(
  titlePanel("SWAT Watershed & Point Source Assessment"),

  sidebarLayout(
    # Convert names to numeric, sort them, then convert back to character for the input

    sidebarPanel(
      # In sidebarPanel
      selectizeInput("select_cach", "Select Outlet Catchment:",
                     choices = NULL, # Leave NULL initially
                     options = list(placeholder = 'Type to search ID...')),
      hr(),
      h4("Statistics for Selection:"),
      verbatimTextOutput("stats_text"),
      downloadButton("download_data", "Export PS Data")
    ),

    mainPanel(
      leafletOutput("map", height = 500),
      hr(),
      tabsetPanel(
        tabPanel("Inflow Data", DTOutput("inflow_table")),
        tabPanel("Point Source Loads", DTOutput("ps_table"))
      )
    )
  )
)

server <- function(input, output, session) {

  # 1. Create a reactive variable to track which ID the tables should show
  # Initialize it as empty
  id_for_tables <- reactiveVal("")

  # 2. Update dropdown choices on start
  updateSelectizeInput(session, "select_cach", choices = choices_list, server = TRUE)

  # 3. When the DROPDOWN changes: Update BOTH the map and the tables
  observeEvent(input$select_cach, {
    req(input$select_cach != "")
    id_for_tables(input$select_cach)
  })

  # 4. When the MAP IS CLICKED: Update ONLY the tables
  observeEvent(input$map_shape_click, {
    req(input$map_shape_click$id)
    # This updates the tables but NOT input$select_cach, so the map won't renew
    id_for_tables(as.character(input$map_shape_click$id))
  })

  # 5. Reactive for Tables (Listens to id_for_tables)
  table_data <- reactive({
    req(id_for_tables() %in% names(river_map))
    river_map[[id_for_tables()]]
  })

  # 6. Reactive for Map (Listens ONLY to the dropdown)
  # This ensures the blue network only changes when the search box is used
  geo_data <- reactive({
    req(input$select_cach)
    req(input$select_cach %in% names(river_map))
    relevant_ids <- river_map[[input$select_cach]]$inflow$cach_id

    list(
      basins   = basins_4326 %>% filter(cach_id %in% relevant_ids),
      segments = segments_4326 %>% filter(id %in% relevant_ids),
      ps       = ps_info %>% filter(cach_id %in% relevant_ids)
    )
  })

  # 7. Static Base Map
  output$map <- renderLeaflet({
    bbox_init <- st_bbox(segments_4326)
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      fitBounds(
        lng1 = as.numeric(bbox_init[1]), lat1 = as.numeric(bbox_init[2]),
        lng2 = as.numeric(bbox_init[3]), lat2 = as.numeric(bbox_init[4])
      ) %>%
      addPolylines(
        data = segments_4326, group = "All Rivers", layerId = ~id,
        weight = 1, color = "#34495e", opacity = 0.4,
        label = ~paste("Segment ID:", id),
        highlightOptions = highlightOptions(weight = 5, color = "red", bringToFront = TRUE),
        labelOptions = labelOptions(sticky = TRUE)
      ) |>
      addLayersControl(

        overlayGroups = c("All Rivers", "Basins", "Rivers", "Point Sources"),

        options = layersControlOptions(collapsed = FALSE)

      )
  })

  # 8. Map Observer (Triggered ONLY by dropdown)
  observe({
    req(geo_data())
    g <- geo_data()
    proxy <- leafletProxy("map") %>%
      clearGroup("Basins") %>% clearGroup("Rivers") %>% clearGroup("Point Sources")

    proxy %>% addPolygons(data = g$basins, group = "Basins", layerId = ~cach_id,
                          weight = 1, color = "grey", fillColor = "blue", fillOpacity = 0.1,
                          popup = ~paste0(
                            "<b>Cach id:</b> ", cach_id, "<br>"))
    proxy %>% addPolylines(data = g$segments, group = "Rivers", layerId = ~id,
                           weight = 3, color = "darkblue", opacity = 0.8)

    # Point Sources
    if (nrow(g$ps) > 0) {
      proxy %>% addCircleMarkers(
        data = g$ps,
        group = "Point Sources",
        radius = 5, color = "red", stroke = FALSE, fillOpacity = 0.8,
        popup = ~paste0(
          "<b>Name:</b> ", name, "<br>",
          "<b>PS CODE:</b> ", ps_code, "<br>",
          "<b>TN Load:</b> ", TN, " kg/y<br>",
          "<b>TP Load:</b> ", TP, " kg/y<br>",
          "<b>Volume:</b> ", volume, " tūkst. m3"
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

  output$ps_table <- renderDT({
    req(table_data())
    datatable(table_data()$ps_load, options = list(pageLength = 5, scrollX = TRUE))
  })
}

shinyApp(ui, server)
