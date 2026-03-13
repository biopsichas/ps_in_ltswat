# Function to read SWAT output files, compute concentrations, and aggregate by month and unit
read_model <- function(f_path, sc_name, subb, setup) {
  dt <- as.data.table(SWATreadR::read_swat(f_path))
  dt[, year_month := sprintf("%04d-%02d", yr, mon)]
  # Compute concentrations
  dt[, tn_conc := ifelse(flo_out == 0, NA_real_,
                         ((orgn_out + nh3_out + no3_out + no2_out)*1e6)/(flo_out*86400*1000))]
  dt[, tp_conc := ifelse(flo_out == 0, NA_real_,
                         ((sedp_out + solp_out)*1e6)/(flo_out*86400*1000))]
  # Add identifiers once
  dt[, `:=`(
    Subbasin = subb,
    Setup_name = setup,
    Scenario = sc_name
  )]
  # Keep only necessary columns
  dt <- dt[, .(unit, year_month, flo_out, tn_conc, tp_conc,
               Subbasin, Setup_name, Scenario)]
  # Ultra-fast aggregation
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE),
           by = .(unit, year_month, Subbasin, Setup_name, Scenario)]

  return(dt)
}

# Function to read reservoirs SWAT files and combine with reservoir.con for GIS ID
read_res <- function(f_path) {
  dt_res <- as.data.table(SWATreadR::read_swat(f_path))|>
    select(name, area_ps, vol_ps, area_es, vol_es)

  as.data.table(SWATreadR::read_swat(sub("hydrology\\.res$", "reservoir.con", f_path))) |>
    select(name, gis_id) |>
    left_join(dt_res, by = "name") |>
    select(gis_id, everything())
}

# Funtion for loading tables
load_table <- function(con, schema, table, exclude_geom = FALSE) {
  # Input validation
  stopifnot(is.character(schema), length(schema) == 1)
  stopifnot(is.character(table), length(table) == 1)
  stopifnot(is.logical(exclude_geom), length(exclude_geom) == 1)

  if (exclude_geom) {
    all_fields <- dbListFields(con, Id(schema = schema, table = table))
    keep_fields <- setdiff(all_fields, c("shape", "geometry"))
    if (length(keep_fields) == 0) {
      stop("All fields excluded, nothing to select.")
    }
    sel_fields <- paste(DBI::dbQuoteIdentifier(con, keep_fields), collapse = ", ")
    message("ⓘ Geometry columns excluded from selection.")
  } else {
    sel_fields <- "*"
  }
  # Quote schema.table properly
  tbl_id <- DBI::dbQuoteIdentifier(con, Id(schema = schema, table = table))

  # Build query
  query <- paste0("SELECT ", sel_fields, " FROM ", tbl_id)

  df <- DBI::dbGetQuery(con, query)
  message("✔ Data is loaded from the database")

  return(df)
}

# Function to get database schema info
get_db_info <- function(con, print_info = TRUE){
  # 1. Get all schemas
  schemas <- dbGetQuery(con, "
  SELECT schema_name
  FROM information_schema.schemata
  WHERE schema_name NOT IN ('pg_catalog', 'information_schema')
")

  # 2. Get all tables
  tables <- dbGetQuery(con, "
  SELECT table_schema, table_name
  FROM information_schema.tables
  WHERE table_type = 'BASE TABLE'
    AND table_schema NOT IN ('pg_catalog', 'information_schema')
")

  # 3. Get all columns and their types
  columns <- dbGetQuery(con, "
  SELECT table_schema, table_name, column_name, data_type
  FROM information_schema.columns
  WHERE table_schema NOT IN ('pg_catalog', 'information_schema')
")

  # 4. Get constraints (PK, FK, UNIQUE, CHECK)
  constraints <- dbGetQuery(con, "
  SELECT
      tc.table_schema,
      tc.table_name,
      kcu.column_name,
      tc.constraint_type
  FROM
      information_schema.table_constraints tc
  LEFT JOIN
      information_schema.key_column_usage kcu
      ON tc.constraint_name = kcu.constraint_name
      AND tc.table_schema = kcu.table_schema
  WHERE
      tc.table_schema NOT IN ('pg_catalog', 'information_schema')
")

  # 5. Get nullability + defaults
  nulls <- dbGetQuery(con, "
  SELECT
      table_schema,
      table_name,
      column_name,
      is_nullable,
      column_default
  FROM
      information_schema.columns
  WHERE
      table_schema NOT IN ('pg_catalog', 'information_schema')
")

  # 6. SAFE row counts
  message("Counting rows (this may take a while)...")

  row_counts <- data.frame(
    table_schema = character(),
    table_name   = character(),
    row_count    = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(tables))) {
    sch <- tables$table_schema[i]
    tbl <- tables$table_name[i]
    full_name <- paste0(sch, ".", tbl)

    # Try a simple, fast estimate first via pg_class (instant, safe)
    rel_tuple_query <- sprintf("
      SELECT schemaname, tablename, n_tup_ins - n_tup_del + n_tup_hot_upd AS estimated_rows
      FROM pg_stat_user_tables
      WHERE schemaname = %s AND tablename = %s
    ", dbQuoteLiteral(con, sch), dbQuoteLiteral(con, tbl))

    est <- tryCatch({
      res <- dbGetQuery(con, rel_tuple_query)
      if (nrow(res) > 0) res$estimated_rows[1] else NA
    }, error = function(e) NA)

    # If no stats or zero, fall back to exact count — but wrap in tryCatch
    if (is.na(est) || est == 0) {
      count_query <- sprintf("SELECT COUNT(*) FROM %s.%s", dbQuoteIdentifier(con, sch), dbQuoteIdentifier(con, tbl))
      exact <- tryCatch({
        dbGetQuery(con, count_query)[1,1]
      }, error = function(e) {
        warning(sprintf("Skipping row count for %s (table inaccessible or corrupted: %s)", full_name, e$message))
        -1  # marker for "failed"
      })
      row_count <- if (is.numeric(exact) && exact >= 0) exact else NA
    } else {
      row_count <- est
    }

    row_counts <- rbind(row_counts, data.frame(
      table_schema = sch,
      table_name   = tbl,
      row_count    = as.numeric(row_count),
      stringsAsFactors = FALSE
    ))

    # Progress feedback
    if (i %% 1 == 0 || i == nrow(tables)) {
      cat(sprintf("\rProcessed %d/%d tables...", i, nrow(tables)))
    }
  }
  cat("\n")

  # 7. Merge everything

  schema_map <- columns %>%
    left_join(nulls, by = c("table_schema", "table_name", "column_name")) %>%
    left_join(constraints, by = c("table_schema", "table_name", "column_name")) %>%
    left_join(row_counts, by = c("table_schema", "table_name")) %>%
    arrange(table_schema, table_name, column_name) %>%
    select(table_schema, table_name, row_count, column_name, data_type,
           is_nullable, column_default, constraint_type)

  # 8. collapse multiple constraint types into one column
  schema_map <- schema_map %>%
    group_by(table_schema, table_name, row_count, column_name,
             data_type, is_nullable, column_default) %>%
    summarise(constraints = paste(unique(constraint_type[!is.na(constraint_type)]), collapse = ", "),
              .groups = "drop")

  # Replace NA row_count with note
  schema_map$row_count[is.na(schema_map$row_count)] <- -1

  # 9. Print summary info if requested
  if (print_info) {
    for (s in unique(schema_map$table_schema)) {
      cat(sprintf("\nSchema: %s\n", s))
      tbls <- schema_map %>% filter(table_schema == s) %>% distinct(table_name, row_count)
      for (i in seq_len(nrow(tbls))) {
        t <- tbls$table_name[i]
        rc <- tbls$row_count[i]
        rc_text <- if (rc == -1) " [inaccessible]" else paste0(" (", scales::comma(rc), " rows)")
        cols <- schema_map %>%
          filter(table_schema == s, table_name == t) %>%
          pull(column_name)
        cat(sprintf("  • %s%s: %s\n", t, rc_text, paste(cols, collapse = ", ")))
      }
    }
  }

  if (nrow(schema_map) == 0) {
    message("No user-defined schemas found.")
  }

  return(schema_map)
}

# Function to write data frame into a table with options for schema creation, overwrite, append, and custom constraints
write_in_table <- function(con, schema, table, df,
                           overwrite = FALSE,
                           append = FALSE,
                           table_constraint = NULL) {
  # Input validation
  stopifnot(is.character(schema), length(schema) == 1)
  stopifnot(is.character(table), length(table) == 1)
  stopifnot(is.data.frame(df))

  # Check schema existence
  schema_exists <- dbGetQuery(con, sprintf("
    SELECT EXISTS(
      SELECT 1
      FROM information_schema.schemata
      WHERE schema_name = '%s'
    ) AS schema_exists;", schema))$schema_exists

  # Create schema if missing
  if (!schema_exists) {
    dbExecute(con, sprintf("CREATE SCHEMA %s;", DBI::dbQuoteIdentifier(con, schema)))
    message(sprintf("✔ Schema '%s' created.", schema))
  }

  # Initialize variable
  colname <- NULL

  # If table_constraint is provided → create table manually
  if (!is.null(table_constraint)) {
    sql <- sprintf(
      "CREATE TABLE IF NOT EXISTS %s.%s (%s);",
      DBI::dbQuoteIdentifier(con, schema),
      DBI::dbQuoteIdentifier(con, table),
      table_constraint
    )

    # Regex to capture: column name + geometry(type, crs)
    matches <- regmatches(
      table_constraint,
      gregexpr("([a-zA-Z0-9_]+)\\s+geometry\\(([^)]*)\\)", table_constraint, perl = TRUE)
    )[[1]]

    # Extract parts
    colname <- gsub(" geometry\\(.*", "", matches)
    # Extract geometry type and CRS
    geom_type <- sub(".*geometry\\(([^,]+),.*", "\\1", matches)
    geom_crs  <- as.integer(sub(".*,(\\s*[0-9]+)\\)$", "\\1", matches))

    if(length(colname) == 0){
      message("ⓘ️ No geometry column found in table_constraint.")
    } else {
      message(paste("ⓘ Geometry definition found:", colname, geom_type, "epsg=", geom_crs, ". It will be used to set geometry type and CRS."))
      if(length(colname) != 1){
        stop("⚠ Multiple geometry columns found in table_constraint. Fix it and try again")
      }
      ## Fix geometry column in df if needed
      # df <- df |>
      #   mutate(!!sym(colname) := st_as_sfc(!!sym(colname), EWKB = TRUE, crs = geom_crs)) |>
      #   st_as_sf() |>
      #   st_make_valid() |>
      #   mutate(!!sym(colname) := st_cast(st_collection_extract(!!sym(colname)), geom_type))
      if(nrow(df) > 0){
        df <- fix_sf_geometry(df, geom_col = colname, geom_crs = geom_crs, geom_type = geom_type)
      }
      # Try to check PostGIS version
      postgis_installed <- TRUE
      tryCatch({
        dbGetQuery(con, "SELECT PostGIS_Version();")
      }, error = function(e) {
        message("⚠ PostGIS not found.")
        postgis_installed <<- FALSE
      })

      # If not installed, create a schema and install PostGIS
      if (!postgis_installed) {
        message("ⓘ  Attempting to create PostGIS extension...")
        # Create public schema if it doesn't exist
        dbExecute(con, "CREATE SCHEMA IF NOT EXISTS public;")
        message("✔ 'public' schema ensured.")

        # Create PostGIS extension in the public schema
        dbExecute(con, "CREATE EXTENSION IF NOT EXISTS postgis WITH SCHEMA public;")
        message("✔ PostGIS extension created in 'public' schema.")
      }

      #  Verify PostGIS is now available
      version <- dbGetQuery(con, "SELECT PostGIS_Version();")
      message(paste("✔ PostGIS version:", version[[1]], "\n"))

    }

    # Create table with constraints
    dbExecute(con, sql)
    message("✔ Table created with custom constraints.")

    if(inherits(df, "sf")){
      rpostgis::pgWriteGeom(
        conn = con,                       # your DBI connection
        name = c(schema, table),          # schema + table as character vector
        data.obj = df,                    # your sf object
        geom = colname,              # name of the geometry column in your table
        overwrite = FALSE,                # append, do NOT overwrite table
        row.names = FALSE                 # do not include R row names
      )

      sql_populate <- sprintf(
        "SELECT Populate_Geometry_Columns('%s.%s'::regclass);",
        DBI::dbQuoteIdentifier(con, schema),
        DBI::dbQuoteIdentifier(con, table)
      )

      dbExecute(con, sql_populate)
      message(sprintf("✔ Geometry column '%s' (%s, SRID %d) added properly.",
                      colname, geom_type, geom_crs))
    } else {
      # Only append, since overwrite would drop constraints
      dbWriteTable(con, Id(schema = schema, table = table), df,
                   append = TRUE, row.names = FALSE)
    }
  } else {
    # Let dbWriteTable handle overwrite/append
    dbWriteTable(con, Id(schema = schema, table = table), df,
                 overwrite = overwrite, append = append, row.names = FALSE)
  }

  message(sprintf("✔ Data written to '%s.%s'.", schema, table))
}

# Function to clean up the database by dropping all user-defined schemas except those specified to skip
vacum_db <- function(con, skip_schemas = c("pg_catalog", "information_schema", "pg_toast")) {
  schemas <- dbGetQuery(con, sprintf("
    SELECT schema_name
    FROM information_schema.schemata
    WHERE schema_name NOT IN (%s)
    ORDER BY schema_name;",
                                     paste0("'", skip_schemas, "'", collapse = ",")
  ))$schema_name

  message("🧹 Cleaning up database...")

  if (length(schemas) == 0) {
    message("ⓘ No schemas to drop.")
    return(invisible(NULL))
  }

  for (s in schemas) {
    tryCatch({
      dbExecute(con, sprintf("DROP SCHEMA IF EXISTS %s CASCADE;", DBI::dbQuoteIdentifier(con, s)))
      message(sprintf("✔ Schema '%s' dropped.", s))
    }, error = function(e) {
      warning(sprintf("⚠ Failed to drop schema '%s': %s", s, e$message))
    })
  }

  message("✔ Database cleanup complete.")
  invisible(NULL)
}

fix_sf_geometry <- function(df, geom_col, geom_crs, geom_type) {
  # ensure column exists
  stopifnot(geom_col %in% names(df))

  # 2 Convert to sf if not already
  if (!inherits(df, "sf")) {

    # If geometry column is ALREADY an sfc object → make sf directly
    if (inherits(df[[geom_col]], "sfc")) {
      df <- sf::st_sf(df, sf_column_name = geom_col, crs = geom_crs)

    } else {
      # Otherwise convert using st_as_sfc (for WKT/WKB)
      df[[geom_col]] <- sf::st_as_sfc(df[[geom_col]], EWKB = TRUE, crs = geom_crs)
      df <- sf::st_sf(df, sf_column_name = geom_col, crs = geom_crs)
    }
  }

  # 3 Ensure geometry column is properly set
  if (is.null(sf::st_geometry(df))) {
    sf::st_geometry(df) <- geom_col
  }

  # 4 Make geometries valid (repair invalid ones)
  df <- sf::st_make_valid(df)

  # 5 Ensure correct geometry type (Polygon → MultiPolygon etc.)
  df[[geom_col]] <- sf::st_cast(df[[geom_col]], toupper(geom_type), warn = FALSE)

  # 6 Return fixed sf object
  return(df)
}

# Function to compare column names and types between two data frames
compare_columns <- function(df1, df2) {
  # Compare names
  name_match <- names(df1) == names(df2)

  # Compare types
  type_match <- sapply(df1, class) == sapply(df2, class)

  # Find mismatches
  mismatch_idx <- which(!name_match | !type_match)

  if(length(mismatch_idx) > 0){
    message("❌ Columns with mismatched names or types:")
    for(i in mismatch_idx){
      message(sprintf("Column %d:", i))
      message("  df1 name/type: ", names(df1)[i], "/", paste(class(df1[[i]]), collapse=", "))
      message("  df2 name/type: ", names(df2)[i], "/", paste(class(df2[[i]]), collapse=", "))
    }
  } else {
    message("✔ All column names and types match.")
  }
}

## Function to read recall input files, process them, and combine with additional data
read_recall_file <- function(f_path, fname, dt3_row) {
  SWATreadR:::read_tbl(file.path(f_path, fname), n_skip = 2) |>
    as_tibble() |>
    mutate(
      date = make_date(yr, mo, day_mo),
      days_in_month = days_in_month(date),
      dplyr::across(all_of(cols_to_scale), \(x) x * days_in_month)
    ) |>
    summarise(
      dplyr::across(all_of(cols_to_scale), sum),
      .by = yr
    ) |>
    bind_cols(dt3_row)
}

# Function to read recall output file for a given subbasin and setup
read_recall_output <- function(subb, setup) {
  base_path <- file.path(setup_path, subb, setup, sc_base,  "recall_yr.txt")

  if (!file.exists(base_path) || length(count.fields(base_path)) < 4) return(NULL)

  match_val <- ps_max$val[ps_max$Subbasin == subb & ps_max$Setup_name == setup]
  ps_count <- if (length(match_val) > 0) match_val[1] else return(NULL)
  SWATreadR::read_swat(base_path) |>
    mutate(id = readr::parse_number(name),
           Subbasin = subb,
           Setup_name = setup) |>
    filter(id <= ps_count) |>
    select(Subbasin, Setup_name, id, yr, all_of(cols_to_scale)) |>
    rename(year = yr)
}

# Reading apportionment excel files
read_apportionment <- function(file, year) {
  read_xlsx(file, sheet = "catchmentloads") |>
    select(CATCHMENTID, ends_with("POINTSOURCE")) |>
    rename(
      gis_id = CATCHMENTID,
      no3  = NO3_POINTSOURCE,
      ntot = NTOTAL_POINTSOURCE,
      solp = PO4_POINTSOURCE,
      ptot = PTOTAL_POINTSOURCE
    ) |>
    mutate(year = year)
}

# Function to check balance and summarize data frame
check_balance <- function(df) {
  # 1. Date Handling: Convert to year if it's a string/date
  if (is.Date(df$year[[1]]) && grepl("^\\d{4}-\\d{2}-\\d{2}$", df$year[[1]])) {
    df <- mutate(df, year = year(as.Date(year)))
  }

  # 2. Aggregation Logic
  df_summary <- df |>
    select(-any_of("id")) |>
    # Aggregate by year and gis_id first, then immediately by year
    group_by(year, gis_id) |>
    summarise(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = "drop_last") |>
    summarise(across(-gis_id, \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
    # Convert to thousands and round
    mutate(across(-c(year), \(x) round(x / 1000, 1))) |>
    # Calculate Mean across all years
    summarise(across(-year, \(x) round(mean(x, na.rm = TRUE), 1)))

  # 3. Column Ordering
  desired_order <- c("flo", "ntot", "ptot", "sed", "orgn", "no3",
                     "no2", "nh3", "solp", "sedp", "cbod")

  df_final <- df_summary |>
    select(starts_with("flo"),
           starts_with("ntot"),
           starts_with("ptot"),
           starts_with("sed"),
           starts_with("orgn"),
           starts_with("no3"),
           starts_with("no2"),
           starts_with("nh3"),
           starts_with("solp"),
           starts_with("sedp"),
           starts_with("cbod"))

  # 4. Feedback
  message("Units: Summarized data in t/year, flo in 1000 m3/year and sed in 1000 t/year.")

  return(df_final)
}
