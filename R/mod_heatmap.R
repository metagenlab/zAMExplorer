mod_heatmap_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(
        title = "Heatmap settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::selectInput(ns("transform"), "Scale", c("Relative abundance (%)" = "relative", "log10(relative abundance + pseudocount)" = "log10p", "CLR" = "clr")),
        shiny::selectInput(ns("annotation"), "Annotate samples by", choices = NULL),
        shiny::numericInput(ns("top_taxa"), "Top taxa to display", value = 30, min = 5, max = 100, step = 5),
        shiny::checkboxInput(ns("cluster_samples"), "Cluster samples", value = TRUE),
        shiny::checkboxInput(ns("cluster_taxa"), "Cluster taxa", value = TRUE),
        shiny::actionButton(ns("generate"), "Generate heatmap", class = "btn-primary")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Publication heatmap", width = 12, status = "primary", solidHeader = TRUE,
        shiny::plotOutput(ns("heatmap"), height = "760px"),
        shiny::fluidRow(
          shiny::column(3, shiny::selectInput(ns("format"), "Format", c("PDF" = "pdf", "PNG (600 dpi)" = "png", "SVG" = "svg"), selected = "pdf")),
          shiny::column(3, shiny::numericInput(ns("width"), "Width (inches)", 10, min = 4, max = 30, step = 0.5)),
          shiny::column(3, shiny::numericInput(ns("height"), "Height (inches)", 8, min = 4, max = 30, step = 0.5)),
          shiny::column(3, shiny::br(), shiny::downloadButton(ns("download"), "Download heatmap", class = "btn-primary"))
        )
      )
    )
  )
}

zamp_heatmap_draw <- function(ht) {
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
}

zamp_export_heatmap <- function(ht, file, format = "pdf", width = 10, height = 8, dpi = 600) {
  format <- tolower(format)
  if (format == "pdf") {
    if (capabilities("cairo")) grDevices::cairo_pdf(file, width = width, height = height) else grDevices::pdf(file, width = width, height = height)
  } else if (format == "svg") {
    grDevices::svg(file, width = width, height = height)
  } else if (format == "png") {
    if (capabilities("cairo")) grDevices::png(file, width = width, height = height, units = "in", res = dpi, type = "cairo") else grDevices::png(file, width = width, height = height, units = "in", res = dpi)
  } else stop("Unsupported heatmap format: ", format)
  on.exit(grDevices::dev.off(), add = TRUE)
  zamp_heatmap_draw(ht)
  invisible(file)
}

mod_heatmap_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps()))
      vars <- zamp_categorical_variables(ps())
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = if ("Genus" %in% ranks) "Genus" else ranks[1])
      shiny::updateSelectInput(session, "annotation", choices = c("None", vars), selected = "None")
    })

    heatmap_obj <- shiny::eventReactive(input$generate, {
      shiny::req(ps(), input$rank)
      x <- zamp_aggregate_rank(ps(), input$rank)
      counts <- zamp_otu_matrix(x, samples_in_rows = TRUE)
      rel <- zamp_relative_matrix(counts)
      means <- colMeans(rel, na.rm = TRUE)
      keep <- names(sort(means, decreasing = TRUE))[seq_len(min(input$top_taxa, length(means)))]
      mat <- rel[, keep, drop = FALSE]
      if (input$transform == "relative") {
        mat <- mat * 100; legend_title <- "Relative\nabundance (%)"
      } else if (input$transform == "log10p") {
        positive <- mat[mat > 0]; pseudo <- if (length(positive)) min(positive) / 2 else 1e-06
        mat <- log10(mat + pseudo); legend_title <- "log10 relative\nabundance"
      } else {
        mat <- zamp_transform_matrix(counts[, keep, drop = FALSE], "clr"); legend_title <- "CLR"
      }
      tax <- zamp_taxon_labels(x, input$rank)
      label_map <- stats::setNames(tax$Taxon, tax$TaxonID)
      colnames(mat) <- make.unique(unname(label_map[colnames(mat)]))
      plot_mat <- t(mat)
      annotation <- NULL
      if (!identical(input$annotation, "None")) {
        meta <- zamp_metadata(x)
        group <- factor(meta[[input$annotation]][match(colnames(plot_mat), meta$Sample)])
        pal <- stats::setNames(zamp_palette(nlevels(group)), levels(group))
        annotation <- ComplexHeatmap::HeatmapAnnotation(Group = group, col = list(Group = pal), annotation_name_gp = grid::gpar(fontface = "bold"))
      }
      col_fun <- grDevices::hcl.colors(101, palette = if (input$transform == "relative") "Blues 3" else "Blue-Red 3")
      ComplexHeatmap::Heatmap(plot_mat, name = legend_title, col = col_fun,
        cluster_rows = isTRUE(input$cluster_taxa), cluster_columns = isTRUE(input$cluster_samples),
        show_column_names = ncol(plot_mat) <= 60, column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 9), top_annotation = annotation, border = TRUE,
        column_title = paste("Taxonomic rank:", input$rank), column_title_gp = grid::gpar(fontface = "bold", fontsize = 12),
        heatmap_legend_param = list(title_gp = grid::gpar(fontface = "bold")))
    }, ignoreInit = TRUE)

    output$heatmap <- shiny::renderPlot({ shiny::req(heatmap_obj()); zamp_heatmap_draw(heatmap_obj()) }, res = 120)
    output$download <- shiny::downloadHandler(
      filename = function() paste0("heatmap_", zamp_safe_stem(input$rank), "_", Sys.Date(), ".", input$format),
      content = function(file) zamp_export_heatmap(heatmap_obj(), file, input$format, input$width, input$height, 600)
    )
    shiny::observeEvent(input$generate, {
      if (isTRUE(auto_save())) {
        shiny::req(heatmap_obj())
        file <- file.path(output_dir(), paste0("heatmap_", zamp_safe_stem(input$rank), "_", base::format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
        zamp_export_heatmap(heatmap_obj(), file, "pdf", 10, 8, 600)
        shiny::showNotification(paste("Saved:", file), type = "message")
      }
    })
  })
}
