mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Read-depth controls", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("group"), "Group read depth by", choices = NULL),
        shiny::numericInput(ns("min_reads"), "Minimum reads per sample", value = 1000, min = 0, step = 100),
        shiny::actionButton(ns("apply_filter"), "Apply read-depth filter", icon = shiny::icon("filter"), class = "btn-primary"),
        shiny::actionButton(ns("reset_filter"), "Reset filter", icon = shiny::icon("rotate-left")),
        shiny::downloadButton(ns("download_filtered"), "Download filtered phyloseq")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Reads across samples", width = 6, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("hist"), height = "420px"),
        zamp_plot_download_ui(ns("hist_download"), 7, 5)
      ),
      shinydashboard::box(title = "Reads across groups", width = 6, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("group_plot"), height = "420px"),
        zamp_plot_download_ui(ns("group_download"), 7, 5)
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Rarefaction", width = 12, status = "primary", solidHeader = TRUE,
        shiny::numericInput(ns("rarefaction_points"), "Number of rarefaction depths", value = 20, min = 5, max = 50),
        shiny::actionButton(ns("make_rarefaction"), "Generate rarefaction curves", class = "btn-primary"),
        plotly::plotlyOutput(ns("rarefaction"), height = "550px"),
        zamp_plot_download_ui(ns("rare_download"), 8, 5.5)
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Sample read counts", width = 12, status = "primary", solidHeader = TRUE,
        DT::DTOutput(ns("reads_table"))
      )
    )
  )
}

mod_qc_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    filtered <- shiny::reactiveVal(NULL)
    current <- shiny::reactive(if (is.null(filtered())) ps() else filtered())

    shiny::observeEvent(ps(), {
      shiny::req(ps())
      shiny::updateSelectInput(session, "group", choices = c("None", zamp_categorical_variables(ps())), selected = "None")
      filtered(NULL)
    })

    shiny::observeEvent(input$apply_filter, {
      shiny::req(ps())
      keep <- phyloseq::sample_sums(ps()) >= input$min_reads
      shiny::validate(shiny::need(any(keep), "No samples remain at this read threshold."))
      x <- phyloseq::prune_samples(keep, ps())
      x <- phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x)
      filtered(x)
      shiny::showNotification(sprintf("QC filter retained %d/%d samples.", phyloseq::nsamples(x), phyloseq::nsamples(ps())), type = "message")
    })
    shiny::observeEvent(input$reset_filter, filtered(NULL))

    reads_df <- shiny::reactive({
      shiny::req(current())
      data.frame(Sample = phyloseq::sample_names(current()), Reads = as.numeric(phyloseq::sample_sums(current())), stringsAsFactors = FALSE)
    })

    hist_plot <- shiny::reactive({
      ggplot2::ggplot(reads_df(), ggplot2::aes(x = Reads)) +
        ggplot2::geom_histogram(bins = 30, fill = "#0072B2", colour = "white", linewidth = 0.25) +
        ggplot2::scale_x_continuous(labels = scales::label_number(big.mark = ",")) +
        ggplot2::labs(title = "Sequencing depth distribution", x = "Reads per sample", y = "Number of samples") +
        zamp_theme()
    })

    group_plot <- shiny::reactive({
      shiny::req(current())
      shiny::validate(shiny::need(!identical(input$group, "None"), "Choose a metadata variable to compare groups."))
      meta <- zamp_metadata(current())
      df <- dplyr::left_join(reads_df(), meta[, c("Sample", input$group), drop = FALSE], by = "Sample")
      names(df)[3] <- "Group"
      df$Group <- factor(df$Group)
      pal <- stats::setNames(zamp_palette(nlevels(df$Group)), levels(df$Group))
      ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Reads, colour = Group)) +
        ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.55) +
        ggplot2::geom_jitter(width = 0.14, size = 2.2, alpha = 0.75) +
        ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
        ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
        ggplot2::labs(title = "Sequencing depth by group", subtitle = input$group, x = NULL, y = "Reads per sample") +
        zamp_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
    })

    rare_plot <- shiny::eventReactive(input$make_rarefaction, {
      shiny::req(current())
      mat <- zamp_otu_matrix(current(), samples_in_rows = TRUE)
      max_depth <- min(rowSums(mat))
      shiny::validate(shiny::need(max_depth >= 2, "Insufficient reads for rarefaction."))
      depths <- unique(round(seq(max(1, max_depth / input$rarefaction_points), max_depth, length.out = input$rarefaction_points)))
      out <- lapply(seq_len(nrow(mat)), function(i) {
        vals <- vapply(depths, function(d) vegan::rarefy(mat[i, ], sample = min(d, sum(mat[i, ]))), numeric(1))
        data.frame(Sample = rownames(mat)[i], Depth = depths, Richness = vals)
      })
      df <- dplyr::bind_rows(out)
      if (!identical(input$group, "None")) {
        meta <- zamp_metadata(current())[, c("Sample", input$group), drop = FALSE]
        names(meta)[2] <- "Group"
        df <- dplyr::left_join(df, meta, by = "Sample")
      } else df$Group <- "All samples"
      df$Group <- factor(df$Group)
      pal <- stats::setNames(zamp_palette(nlevels(df$Group)), levels(df$Group))
      ggplot2::ggplot(df, ggplot2::aes(x = Depth, y = Richness, group = Sample, colour = Group)) +
        ggplot2::geom_line(alpha = 0.7, linewidth = 0.7) +
        ggplot2::scale_colour_manual(values = pal) +
        ggplot2::labs(title = "Rarefaction curves", x = "Subsampled reads", y = "Expected richness", colour = if (identical(input$group, "None")) NULL else input$group) +
        zamp_theme()
    }, ignoreInit = TRUE)

    output$hist <- plotly::renderPlotly(zamp_plotly(hist_plot()))
    output$group_plot <- plotly::renderPlotly(zamp_plotly(group_plot()))
    output$rarefaction <- plotly::renderPlotly({ shiny::req(rare_plot()); zamp_plotly(rare_plot()) })
    output$reads_table <- DT::renderDT(DT::datatable(reads_df(), options = list(pageLength = 12), rownames = FALSE))

    zamp_plot_download_server("hist_download", hist_plot, function() paste0("reads_distribution_", Sys.Date()))
    zamp_plot_download_server("group_download", group_plot, function() paste0("reads_by_group_", Sys.Date()))
    zamp_plot_download_server("rare_download", rare_plot, function() paste0("rarefaction_", Sys.Date()))

    output$download_filtered <- shiny::downloadHandler(
      filename = function() paste0("filtered_phyloseq_", Sys.Date(), ".rds"),
      content = function(file) saveRDS(current(), file)
    )

    shiny::observeEvent(input$make_rarefaction, {
      if (isTRUE(auto_save())) {
        shiny::req(rare_plot())
        file <- zamp_save_plot_to_dir(rare_plot(), output_dir(), "rarefaction", "pdf", 8, 5.5)
        shiny::showNotification(paste("Saved:", file), type = "message")
      }
    })

    current
  })
}
