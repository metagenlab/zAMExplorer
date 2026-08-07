# Publication-quality plotting and export helpers -----------------------------

zamp_theme <- function(base_size = 12, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, colour = "grey30", margin = ggplot2::margin(b = 8)),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

zamp_palette <- function(n) {
  stopifnot(n >= 0)
  if (n == 0) return(character())
  okabe_ito <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000"
  )
  if (n <= length(okabe_ito)) return(okabe_ito[seq_len(n)])
  c(okabe_ito, grDevices::hcl.colors(n - length(okabe_ito), palette = "Dynamic"))
}

zamp_named_palette <- function(values, other_label = "Other") {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  main <- setdiff(values, other_label)
  pal <- stats::setNames(zamp_palette(length(main)), main)
  if (other_label %in% values) pal <- c(pal, stats::setNames("#BDBDBD", other_label))
  pal
}

zamp_safe_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

zamp_export_plot <- function(plot, file, format = "pdf", width = 7, height = 5, dpi = 600) {
  format <- tolower(format)
  if (!format %in% c("pdf", "png", "svg")) {
    stop("Unsupported plot format: ", format)
  }

  device <- switch(
    format,
    png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png",
    svg = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg",
    pdf = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  )

  ggplot2::ggsave(
    filename = file,
    plot = plot,
    device = device,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  invisible(file)
}

zamp_save_plot_to_dir <- function(plot, output_dir, stem, format = "pdf", width = 7, height = 5, dpi = 600) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(
    output_dir,
    paste0(zamp_safe_stem(stem), "_", base::format(Sys.time(), "%Y%m%d_%H%M%S"), ".", format)
  )
  zamp_export_plot(plot, file, format = format, width = width, height = height, dpi = dpi)
  file
}

zamp_save_table_to_dir <- function(x, output_dir, stem, ext = "tsv") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(
    output_dir,
    paste0(zamp_safe_stem(stem), "_", base::format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
  )
  if (ext == "csv") {
    utils::write.csv(x, file, row.names = FALSE)
  } else {
    utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  file
}

zamp_plot_download_ui <- function(id, default_width = 7, default_height = 5, default_format = "pdf") {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(3, shiny::selectInput(ns("format"), "Format", c("PDF" = "pdf", "PNG (600 dpi)" = "png", "SVG" = "svg"), selected = default_format)),
    shiny::column(3, shiny::numericInput(ns("width"), "Width (inches)", default_width, min = 3, max = 30, step = 0.5)),
    shiny::column(3, shiny::numericInput(ns("height"), "Height (inches)", default_height, min = 3, max = 30, step = 0.5)),
    shiny::column(3, shiny::br(), shiny::downloadButton(ns("download"), "Download publication plot", class = "btn-primary"))
  )
}

zamp_plot_download_server <- function(id, plot_reactive, filename_stem) {
  shiny::moduleServer(id, function(input, output, session) {
    output$download <- shiny::downloadHandler(
      filename = function() paste0(zamp_safe_stem(filename_stem()), ".", input$format),
      content = function(file) {
        shiny::req(plot_reactive())
        zamp_export_plot(
          plot_reactive(), file,
          format = input$format,
          width = input$width,
          height = input$height,
          dpi = 600
        )
      }
    )
  })
}

zamp_plotly <- function(plot, tooltip = "all") {
  plotly::ggplotly(plot, tooltip = tooltip, dynamicTicks = FALSE) |>
    plotly::config(
      displaylogo = FALSE,
      scrollZoom = TRUE,
      modeBarButtonsToRemove = c("lasso2d", "select2d")
    )
}
