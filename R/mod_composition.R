mod_composition_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Composition settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::numericInput(ns("detection"), "Minimum abundance in any sample (%)", 0.1, min = 0, max = 100, step = 0.05),
        shiny::numericInput(ns("prevalence"), "Minimum prevalence (%)", 10, min = 0, max = 100, step = 1),
        shiny::selectInput(ns("order_by"), "Order samples by metadata", choices = NULL),
        shiny::selectInput(ns("facet_by"), "Facet by metadata", choices = NULL),
        shiny::actionButton(ns("generate"), "Generate composition plot", class = "btn-primary")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Relative abundance", width = 12, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("plot"), height = "720px"),
        zamp_plot_download_ui(ns("download"), 10, 6)
      )
    )
  )
}

mod_composition_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps()))
      default_rank <- if ("Genus" %in% ranks) "Genus" else ranks[min(2, length(ranks))]
      meta <- zamp_all_metadata_variables(ps())
      group_vars <- zamp_categorical_variables(ps())
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = default_rank)
      shiny::updateSelectInput(session, "order_by", choices = c("None", meta), selected = "None")
      shiny::updateSelectInput(session, "facet_by", choices = c("None", group_vars), selected = "None")
    })

    composition_data <- shiny::eventReactive(input$generate, {
      shiny::req(ps(), input$rank)
      x <- zamp_aggregate_rank(ps(), input$rank)
      mat <- zamp_otu_matrix(x, samples_in_rows = TRUE)
      rel <- zamp_relative_matrix(mat)
      labels <- zamp_rank_labels(x, input$rank)
      colnames(rel) <- labels

      prevalence <- colMeans(rel > (input$detection / 100)) * 100
      max_abundance <- apply(rel, 2, max) * 100
      keep <- prevalence >= input$prevalence & max_abundance >= input$detection
      if (!any(keep)) keep[which.max(colMeans(rel))] <- TRUE

      kept <- rel[, keep, drop = FALSE]
      other <- if (any(!keep)) rowSums(rel[, !keep, drop = FALSE]) else rep(0, nrow(rel))
      if (any(other > 0)) kept <- cbind(kept, Other = other)

      df <- as.data.frame(kept, check.names = FALSE)
      df$Sample <- rownames(df)
      long <- tidyr::pivot_longer(df, -Sample, names_to = "Taxon", values_to = "Abundance")
      meta <- zamp_metadata(x)
      long <- dplyr::left_join(long, meta, by = "Sample")

      tax_order <- names(sort(tapply(long$Abundance, long$Taxon, mean), decreasing = TRUE))
      if ("Other" %in% tax_order) tax_order <- c(setdiff(tax_order, "Other"), "Other")
      long$Taxon <- factor(long$Taxon, levels = tax_order)

      sample_order <- unique(long$Sample)
      if (!identical(input$order_by, "None") && input$order_by %in% names(long)) {
        ord_df <- unique(long[, c("Sample", input$order_by), drop = FALSE])
        ord_df <- ord_df[order(ord_df[[input$order_by]], ord_df$Sample, na.last = TRUE), , drop = FALSE]
        sample_order <- ord_df$Sample
      }
      long$Sample <- factor(long$Sample, levels = sample_order)
      long
    }, ignoreInit = TRUE)

    comp_plot <- shiny::reactive({
      df <- composition_data()
      pal <- zamp_named_palette(levels(df$Taxon))
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Abundance, fill = Taxon)) +
        ggplot2::geom_col(width = 0.92, colour = NA) +
        ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
        ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.01))) +
        ggplot2::labs(title = "Microbial community composition", subtitle = paste("Taxonomic rank:", input$rank), x = NULL, y = "Relative abundance", fill = input$rank) +
        zamp_theme(11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), legend.position = "right")
      if (!identical(input$facet_by, "None") && input$facet_by %in% names(df)) {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste("~", input$facet_by)), scales = "free_x", space = "free_x")
      }
      p
    })

    output$plot <- plotly::renderPlotly({ shiny::req(comp_plot()); zamp_plotly(comp_plot()) })
    zamp_plot_download_server("download", comp_plot, function() paste0("composition_", input$rank, "_", Sys.Date()))

    shiny::observeEvent(input$generate, {
      if (isTRUE(auto_save())) {
        shiny::req(comp_plot())
        file <- zamp_save_plot_to_dir(comp_plot(), output_dir(), paste0("composition_", input$rank), "pdf", 10, 6)
        shiny::showNotification(paste("Saved:", file), type = "message")
      }
    })
  })
}
