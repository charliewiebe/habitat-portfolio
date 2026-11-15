# Kit fox habitat & reintroduction outlook: Shiny + leaflet app
# Run from the project root after the pipeline: shiny::runApp("shiny")
ROOT <- normalizePath("..")
local({ owd <- setwd(ROOT); on.exit(setwd(owd)); source("pipeline/00_config.R", local = globalenv()) })
source(file.path(ROOT, "pipeline", "functions.R"))
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(leaflet); library(plotly); library(DT); library(sf)
})

OUT <- file.path(ROOT, CFG$out_dir)
if (!all(file.exists(file.path(OUT, c("window_static.tif", "timeseries.csv", "rangewide_index.tif")))))
  stop("Pipeline outputs not found. Run Rscript run_all.R from the project root first.")

SET_LABELS <- c("All records (effort-corrected)" = "all", "Region-balanced" = "balanced",
                "Exclude San Joaquin" = "no_sjv", "All records, no effort correction" = "all_nobias")
SSP_LABELS <- setNames(paste0("ssp", CFG$ssps),
                       c("245" = "SSP2-4.5", "126" = "SSP1-2.6", "370" = "SSP3-7.0", "585" = "SSP5-8.5")[CFG$ssps])
PERIOD_LABELS <- c("Baseline (1970-2000)" = "baseline", setNames(CFG$periods, CFG$periods))
LAYER_LABELS <- c("Habitat index" = "index", "Suitable habitat & dispersal networks" = "networks",
                  "GCM agreement on suitability" = "agree", "Extrapolation (share of GCMs, MESS < 0)" = "extrap",
                  "Cropland + built fraction" = "unavail", "Population density (effort covariate)" = "pop")
METRICS <- c("Largest network carrying capacity (foxes)" = "largest_network_K",
             "Suitable habitat in Oregon (km²)" = "suitable_km2",
             "Core habitat in Oregon (km²)" = "core_km2",
             "Suitable habitat in Harney + Malheur (km²)" = "hist_range_km2",
             "Networks at or above MVP" = "n_viable",
             "Potential habitat before land-use loss (km²)" = "potential_km2")

static <- rast(file.path(OUT, "window_static.tif"))
idx_stacks  <- lapply(setNames(CFG$model_sets, CFG$model_sets),
                      function(s) rast(file.path(OUT, sprintf("window_index_%s.tif", s))))
mess_stacks <- lapply(setNames(CFG$model_sets, CFG$model_sets),
                      function(s) rast(file.path(OUT, sprintf("window_mess_%s.tif", s))))
rangewide <- rast(file.path(OUT, "rangewide_index.tif"))
effort    <- rast(file.path(OUT, "rangewide_effort.tif"))
ts0  <- fread(file.path(OUT, "timeseries.csv"))
cv   <- fread(file.path(OUT, "cv_summary.csv"))
imp  <- fread(file.path(OUT, "importance.csv"))
rc   <- fread(file.path(OUT, "response_curves.csv"))
orc  <- fread(file.path(OUT, "oregon_check.csv"))
clog <- fread(file.path(ROOT, CFG$data_dir, "occ_cleaning_log.csv"))
md   <- readRDS(file.path(ROOT, CFG$data_dir, "model_data.rds"))
occ  <- md$occ
gpkg <- file.path(ROOT, CFG$data_dir, "env", "boundaries.gpkg")
adm_sf    <- st_as_sf(vect(gpkg, layer = "admin1"))
oregon_sf <- adm_sf[adm_sf$NAME_1 == "Oregon", ]
hist_sf   <- st_as_sf(vect(gpkg, layer = "hist_range"))
or_obs    <- unique(orc[, .(id, observed_on, lon, lat, obscured, url)])

pal_idx <- colorNumeric("viridis", c(0, 1), na.color = "transparent")
fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")

# ---------------------------------------------------------------------------------------------
ui <- page_navbar(
  title = "Kit fox habitat outlook, Oregon",
  theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_collection("system-ui", "sans-serif")),
  nav_panel("Oregon outlook",
    layout_sidebar(
      sidebar = sidebar(width = 320,
        selectInput("set", "Training set", SET_LABELS, selected = "no_sjv"),
        selectInput("ssp", "Emissions scenario", SSP_LABELS),
        selectInput("period", "Period", PERIOD_LABELS),
        selectInput("layer", "Map layer", LAYER_LABELS, selected = "networks"),
        checkboxInput("show_obs", "Show iNaturalist records", TRUE),
        accordion(open = FALSE,
          accordion_panel("Viability assumptions",
            sliderInput("suitable_q", "Suitable if habitat index ≥", 0.05, 0.5, CFG$suitable_q, 0.05),
            sliderInput("core_q", "Core if habitat index ≥", 0.2, 0.9, CFG$core_q, 0.05),
            sliderInput("density", "Density in suitable habitat (foxes/km²)", 0.02, 0.3, CFG$density_km2, 0.01),
            sliderInput("mvp", "Minimum viable population (foxes)", 25, 500, CFG$mvp_n, 25),
            sliderInput("dispersal", "Dispersal distance linking patches (km)", 0, 50, CFG$dispersal_km, 5),
            sliderInput("unavail", "Unavailable if cropland + built ≥", 0.1, 1, CFG$unavailable_frac, 0.1),
            checkboxInput("excl_extrap", "Exclude extrapolated cells (MESS < 0)", CFG$exclude_extrapolation)
          )
        ),
        actionButton("recompute", "Recompute outlook with these assumptions", class = "btn-primary"),
        uiOutput("ts_note")
      ),
      layout_columns(fill = FALSE, col_widths = c(3, 3, 3, 3),
        value_box("Suitable habitat in Oregon", textOutput("vb_suit"), p(textOutput("vb_suit_sub"))),
        value_box("Largest connected network", textOutput("vb_k"), p(textOutput("vb_k_sub"))),
        value_box("Networks ≥ MVP", textOutput("vb_viable"), p(textOutput("vb_viable_sub"))),
        value_box("Reintroduction verdict", textOutput("vb_verdict"), p(textOutput("vb_verdict_sub")))
      ),
      card(full_screen = TRUE, card_header(textOutput("map_title")), leafletOutput("map", height = 520)),
      card(full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Outlook over time (line = GCM median, band = GCM range)"),
          div(style = "width: 340px; margin-bottom: -1rem;",
              selectInput("metric", NULL, METRICS, width = "100%")))),
        plotlyOutput("ts_plot", height = 360))
    )
  ),
  nav_panel("Range-wide model",
    layout_sidebar(
      sidebar = sidebar(width = 300,
        selectInput("rw_set", "Training set", SET_LABELS),
        radioButtons("rw_layer", "Layer", c("Habitat index (effort held constant)" = "index",
                                            "Observer-effort multiplier (log10)" = "effort")),
        p(class = "small text-muted",
          "The effort multiplier is how many times more records the fitted model expects in a cell than",
          "it would at the reference effort level (median population density and travel time of the background).",
          "It comes from the all-records model and is removed before projection.")
      ),
      card(full_screen = TRUE, leafletOutput("rw_map", height = 640))
    )
  ),
  nav_panel("Diagnostics",
    layout_columns(col_widths = c(6, 6, 6, 6, 12, 6, 6),
      card(card_header("Spatial block cross-validation (250 km blocks)"), DTOutput("cv_tbl")),
      card(card_header("Permutation importance (drop in ensemble AUC)"),
           selectInput("imp_set", NULL, SET_LABELS, selected = "no_sjv"), plotlyOutput("imp_plot", height = 300)),
      card(card_header("Response curves (relative intensity; other predictors at presence medians)"),
           layout_columns(selectInput("rc_set", NULL, SET_LABELS, selected = "no_sjv"),
                          selectInput("rc_var", NULL, setNames(md$vars, var_label(md$vars)))),
           plotlyOutput("rc_plot", height = 300)),
      card(card_header("External check: Oregon iNaturalist records (never used in training)"), DTOutput("or_tbl")),
      card(card_header("iNaturalist records by year (cleaned, all states)"), plotlyOutput("yr_plot", height = 250)),
      card(card_header("Record cleaning"), tableOutput("clog_tbl")),
      card(card_header("Presence cells by state (training)"), tableOutput("state_tbl"))
    )
  ),
  nav_panel("Methods", card(card_body(class = "lh-lg", HTML(paste(readLines(file.path(ROOT, "shiny", "methods.html")), collapse = "\n")))))
)

# ---------------------------------------------------------------------------------------------
server <- function(input, output, session) {

  P <- reactive({
    validate(need(input$core_q >= input$suitable_q, "Core threshold must be at least the suitable threshold."))
    modifyList(CFG, list(suitable_q = input$suitable_q, core_q = input$core_q, density_km2 = input$density,
                         mvp_n = input$mvp, dispersal_km = input$dispersal, unavailable_frac = input$unavail,
                         exclude_extrapolation = input$excl_extrap))
  })
  pkey <- reactive(paste(unlist(P()[c("suitable_q", "core_q", "density_km2", "mvp_n", "dispersal_km",
                                      "unavailable_frac", "exclude_extrapolation")]), collapse = "|"))
  default_key <- paste(unlist(CFG[c("suitable_q", "core_q", "density_km2", "mvp_n", "dispersal_km",
                                    "unavailable_frac", "exclude_extrapolation")]), collapse = "|")

  surface <- reactive({
    st <- idx_stacks[[input$set]]
    L <- if (input$period == "baseline") "baseline" else
      grep(paste0("^", input$ssp, "__", input$period, "__"), names(st), value = TRUE)
    validate(need(length(L) > 0, "No layers for this scenario."))
    idx <- st[[L]]; ms <- mess_stacks[[input$set]][[L]]
    list(mean = if (length(L) > 1) mean(idx) else idx,
         agree = mean(idx >= input$suitable_q),
         extrap = mean(ms < 0), n = length(L))
  })
  via <- reactive({
    s <- surface()
    viability(s$mean, static[["unavail"]], s$extrap < 0.5, static[["oregon"]],
              static[["hist_range"]], static[["area_km2"]], P())
  })

  # ---- value boxes ----
  output$vb_suit <- renderText(paste(fmt(via()$metrics$suitable_km2), "km²"))
  output$vb_suit_sub <- renderText(sprintf("Core %s km² · Harney+Malheur %s km²",
                                           fmt(via()$metrics$core_km2), fmt(via()$metrics$hist_range_km2)))
  output$vb_k <- renderText(paste("≈", fmt(via()$metrics$largest_network_K), "foxes"))
  output$vb_k_sub <- renderText(sprintf("%s km² at %.2f foxes/km²", fmt(via()$metrics$largest_network_km2), input$density))
  output$vb_viable <- renderText(fmt(via()$metrics$n_viable))
  output$vb_viable_sub <- renderText(sprintf("of %d networks (MVP %d, dispersal %d km)",
                                             via()$metrics$n_networks, input$mvp, input$dispersal))
  output$vb_verdict <- renderText(if (via()$metrics$n_viable > 0) "Enough habitat" else "Not enough")
  output$vb_verdict_sub <- renderText({
    m <- via()$metrics
    if (m$n_viable > 0) sprintf("%d network(s) hold ≥ %d foxes", m$n_viable, input$mvp)
    else sprintf("largest holds ≈ %s, needs %d", fmt(m$largest_network_K), input$mvp)
  })
  output$map_title <- renderText({
    per <- names(PERIOD_LABELS)[PERIOD_LABELS == input$period]
    sc <- if (input$period == "baseline") "" else paste0(names(SSP_LABELS)[SSP_LABELS == input$ssp],
                                                        sprintf(", mean of %d GCMs", surface()$n))
    paste0(names(LAYER_LABELS)[LAYER_LABELS == input$layer], " · ", per, if (nzchar(sc)) paste0(" · ", sc) else "",
           " · Oregon totals only")
  })

  # ---- Oregon map ----
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addProviderTiles(providers$CartoDB.Positron, group = "Light") |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Imagery") |>
      setView(-119.3, 43.6, 6) |>
      addPolylines(data = st_cast(oregon_sf, "MULTILINESTRING"), color = "#222", weight = 2, group = "outline") |>
      addPolylines(data = st_cast(hist_sf, "MULTILINESTRING"), color = "#6a3d9a", weight = 2, dashArray = "6 4",
                   group = "outline", label = ~paste(NAME_2, "County (historical range)")) |>
      addLayersControl(baseGroups = c("Light", "Imagery"), position = "topright")
  })

  observe({
    s <- surface(); v <- via(); lay <- input$layer
    px <- leafletProxy("map") |> clearGroup("surface") |> clearGroup("networks") |> removeControl("legend")
    if (lay == "networks") {
      if (nrow(v$tab)) {
        net <- v$net; names(net) <- "net"
        vr <- subst(net, v$tab$net, ifelse(v$tab$viable, 2, 1))
        px <- px |> addRasterImage(vr, colors = colorFactor(c("#e08214", "#1b7837"), 1:2, na.color = "transparent"),
                                   opacity = 0.75, group = "surface", method = "ngb")
        polys <- st_as_sf(as.polygons(net)); polys <- merge(polys, v$tab, by = "net")
        px <- px |> addPolygons(data = polys, fill = FALSE, weight = 1, color = "#333", group = "networks",
          popup = ~sprintf("<b>Network %d</b><br>%s km² suitable<br>K ≈ %s foxes<br>%s", net, fmt(km2), fmt(K),
                           ifelse(viable, "≥ MVP", "below MVP")))
      }
      px |> addLegend("bottomright", colors = c("#1b7837", "#e08214"), layerId = "legend",
                      labels = c("Network ≥ MVP", "Network below MVP"), title = "Suitable, available habitat")
    } else {
      spec <- switch(lay,
        index   = list(r = s$mean, pal = pal_idx, title = "Habitat index"),
        agree   = list(r = s$agree, pal = colorNumeric("PuBuGn", c(0, 1), na.color = "transparent"), title = "Share of GCMs suitable"),
        extrap  = list(r = s$extrap, pal = colorNumeric("Reds", c(0, 1), na.color = "transparent"), title = "Share of GCMs extrapolating"),
        unavail = list(r = static[["unavail"]], pal = colorNumeric("YlOrBr", c(0, 1), na.color = "transparent"), title = "Cropland + built"),
        pop     = list(r = static[["log_pop"]], pal = colorNumeric("magma", range(values(static[["log_pop"]]), na.rm = TRUE),
                                                                    na.color = "transparent", reverse = TRUE), title = "log(1 + people/km²)"))
      px |> addRasterImage(spec$r, colors = spec$pal, opacity = 0.75, group = "surface", method = "bilinear") |>
        addLegend("bottomright", pal = spec$pal, values = values(spec$r), title = spec$title, layerId = "legend")
    }
  })

  observe({
    px <- leafletProxy("map") |> clearGroup("obs")
    if (!isTRUE(input$show_obs)) return()
    near <- occ[in_or == FALSE & lon > CFG$window[1] & lon < CFG$window[2] & lat > CFG$window[3] & lat < CFG$window[4]]
    precise <- or_obs[obscured == FALSE]; obsc <- or_obs[obscured == TRUE]
    px <- px |> addCircleMarkers(data = near, ~lon, ~lat, radius = 3, stroke = FALSE, fillColor = "#555",
                                 fillOpacity = 0.8, group = "obs", label = ~paste("Training record", observed_on))
    if (nrow(obsc)) {
      b <- obsc[, .(x0 = floor(lon / 0.2) * 0.2, y0 = floor(lat / 0.2) * 0.2, observed_on)]
      px <- px |> addRectangles(data = b, ~x0, ~y0, ~x0 + 0.2, ~y0 + 0.2, color = "#d7191c", weight = 2,
                                fill = FALSE, group = "obs", label = ~paste("Oregon record (obscured to 0.2°)", observed_on))
    }
    if (nrow(precise))
      px |> addCircleMarkers(data = precise, ~lon, ~lat, radius = 6, color = "#fff", weight = 2, fillColor = "#d7191c",
                             fillOpacity = 1, group = "obs", label = ~paste("Oregon record", observed_on))
  })

  # ---- outlook time series ----
  cache <- reactiveValues()
  for (s in CFG$model_sets) cache[[s]] <- list(key = default_key, data = ts0[set == s])
  observeEvent(input$recompute, {
    withProgress(message = "Recomputing outlook", value = 0.3, {
      d <- outlook_series(idx_stacks[[input$set]], mess_stacks[[input$set]], static, P())
      d[, `:=`(set = input$set, year = CFG$period_mid[ifelse(ssp == "baseline", "baseline", period)])]
      cache[[input$set]] <- list(key = pkey(), data = d)
    })
  })
  ts_cur <- reactive({ req(cache[[input$set]]); cache[[input$set]] })
  output$ts_note <- renderUI({
    if (ts_cur()$key != pkey())
      p(class = "small text-warning mt-2", "Map reflects current assumptions; the time series still uses the previous ones.")
  })
  output$ts_plot <- renderPlotly({
    d <- ts_cur()$data; m <- input$metric
    base <- d[ssp == "baseline"][[m]]
    p <- plot_ly()
    cols <- c(ssp126 = "#1a9850", ssp245 = "#4575b4", ssp370 = "#fdae61", ssp585 = "#d73027")
    for (s in unique(d[ssp != "baseline", ssp])) {
      g <- d[ssp == s, .(med = median(get(m)), lo = min(get(m)), hi = max(get(m))), by = year][order(year)]
      g <- rbind(data.table(year = CFG$period_mid[["baseline"]], med = base, lo = base, hi = base), g)
      lab <- names(SSP_LABELS)[SSP_LABELS == s]
      p <- p |> add_ribbons(data = g, x = ~year, ymin = ~lo, ymax = ~hi, fillcolor = cols[[s]], opacity = 0.18,
                            line = list(width = 0), name = paste(lab, "GCM range"), showlegend = FALSE) |>
        add_lines(data = g, x = ~year, y = ~med, line = list(color = cols[[s]], width = 2.5), name = lab)
    }
    shapes <- if (m == "largest_network_K") list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                y0 = P()$mvp_n, y1 = P()$mvp_n, line = list(dash = "dash", color = "#666"))) else list()
    p |> layout(xaxis = list(title = "", tickvals = unname(CFG$period_mid),
                             ticktext = c("1970-2000", CFG$periods)),
                yaxis = list(title = names(METRICS)[METRICS == m], rangemode = "tozero"),
                shapes = shapes, legend = list(orientation = "h", y = -0.15), hovermode = "x unified")
  })

  # ---- range-wide map ----
  output$rw_map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addProviderTiles(providers$CartoDB.Positron) |> setView(-111, 35, 5) |>
      addPolylines(data = st_cast(oregon_sf, "MULTILINESTRING"), color = "#222", weight = 2)
  })
  observe({
    px <- leafletProxy("rw_map") |> clearGroup("surface") |> clearGroup("obs") |> removeControl("legend")
    if (input$rw_layer == "index") {
      px <- px |> addRasterImage(rangewide[[input$rw_set]], colors = pal_idx, opacity = 0.7, group = "surface") |>
        addLegend("bottomright", pal = pal_idx, values = c(0, 1), title = "Habitat index", layerId = "legend")
    } else {
      rng <- range(values(effort), na.rm = TRUE)
      pal <- colorNumeric("RdBu", c(-max(abs(rng)), max(abs(rng))), reverse = TRUE, na.color = "transparent")
      px <- px |> addRasterImage(effort, colors = pal, opacity = 0.7, group = "surface") |>
        addLegend("bottomright", pal = pal, values = rng, title = "log10 effort multiplier", layerId = "legend")
    }
    used <- if (input$rw_set == "no_sjv") !occ$in_sjv & !occ$in_or else !occ$in_or
    px |> addCircleMarkers(data = occ[used], ~lon, ~lat, radius = 3, stroke = FALSE, fillColor = "#111",
                           fillOpacity = 0.7, group = "obs", label = ~paste(state, observed_on)) |>
      addCircleMarkers(data = occ[!used], ~lon, ~lat, radius = 4, color = "#d7191c", weight = 1.5, fill = FALSE,
                       group = "obs", label = ~ifelse(in_or, "Oregon record (held out)", "Excluded from this set"))
  })

  # ---- diagnostics ----
  output$cv_tbl <- renderDT({
    d <- cv[, .(Set = names(SET_LABELS)[match(set, SET_LABELS)], Algorithm = algo, AUC = round(auc, 3),
                `AUC sd` = round(auc_sd, 3), Boyce = round(boyce, 3), `Boyce sd` = round(boyce_sd, 3),
                `Held-out presences ≥ P10` = round(test_pres_above_p10, 3))]
    datatable(d, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })
  output$imp_plot <- renderPlotly({
    d <- imp[set == input$imp_set][order(share)][, label := var_label(variable)]
    plot_ly(d, x = ~share, y = ~factor(label, levels = label), type = "bar", marker = list(color = "#2c7fb8")) |>
      layout(xaxis = list(title = "Share of AUC drop", tickformat = ".0%"), yaxis = list(title = ""))
  })
  output$rc_plot <- renderPlotly({
    d <- rc[set == input$rc_set & variable == input$rc_var]
    plot_ly(d, x = ~x, y = ~relative, color = ~algo, colors = c("#1b9e77", "#7570b3"), type = "scatter", mode = "lines") |>
      layout(xaxis = list(title = var_label(input$rc_var)), yaxis = list(title = "Relative intensity", range = c(0, 1.05)))
  })
  output$or_tbl <- renderDT({
    d <- dcast(orc, id + observed_on + obscured ~ set, value.var = "index_mean")
    num <- intersect(CFG$model_sets, names(d)); d[, (num) := lapply(.SD, round, 2), .SDcols = num]
    datatable(d, rownames = FALSE, options = list(dom = "t"),
              caption = "Habitat index at each record (obscured records: mean over their 0.2° cell). ≥ 0.10 = suitable.")
  })
  output$yr_plot <- renderPlotly({
    d <- occ[, .N, by = year][order(year)]
    plot_ly(d, x = ~year, y = ~N, type = "bar", marker = list(color = "#636363")) |>
      layout(xaxis = list(title = ""), yaxis = list(title = "Records"))
  })
  output$clog_tbl <- renderTable(clog, striped = TRUE)
  output$state_tbl <- renderTable(md$n_state[, .(State = state, `Presence cells` = N)], striped = TRUE)
}

shinyApp(ui, server)
