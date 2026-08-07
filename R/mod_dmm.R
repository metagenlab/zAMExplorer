mod_dmm_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Dirichlet multinomial mixture settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::numericInput(ns("detection"), "Detection threshold (%)", 0.1, min = 0, max = 100, step = 0.05),
        shiny::numericInput(ns("prevalence"), "Prevalence threshold (%)", 10, min = 0, max = 100, step = 1),
        shiny::numericInput(ns("max_k"), "Maximum number of components", 6, min = 2, max = 12),
        shiny::actionButton(ns("run"), "Fit DMM models", class = "btn-primary")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Model selection", width = 6, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("fit_plot"), height = "450px"), zamp_plot_download_ui(ns("fit_download"), 6.5, 5)),
      shinydashboard::box(title = "Community assignments", width = 6, status = "primary", solidHeader = TRUE,
        shiny::verbatimTextOutput(ns("best_model")), DT::DTOutput(ns("assignments")),
        shiny::downloadButton(ns("assignments_download"), "Download assignments"),
        shiny::downloadButton(ns("phyloseq_download"), "Download phyloseq with DMM cluster"))
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Community-type drivers", width = 12, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("drivers"), height = "700px"), zamp_plot_download_ui(ns("drivers_download"), 10, 7))
    )
  )
}

mod_dmm_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps()); ranks <- c("ASV", phyloseq::rank_names(ps()))
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = if ("Genus" %in% ranks) "Genus" else ranks[1])
    })
    fitted_models <- shiny::eventReactive(input$run, {
      shiny::req(ps(), input$rank)
      x <- zamp_aggregate_rank(ps(), input$rank); counts <- zamp_otu_matrix(x, samples_in_rows = TRUE); rel <- zamp_relative_matrix(counts)
      keep <- colMeans(rel > input$detection / 100) * 100 >= input$prevalence
      shiny::validate(shiny::need(any(keep), "No taxa remain after detection/prevalence filtering."))
      counts <- counts[, keep, drop = FALSE]; counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      shiny::validate(shiny::need(nrow(counts) >= 4, "At least four non-empty samples are required for DMM."))
      max_k <- min(input$max_k, nrow(counts) - 1)
      fits <- lapply(seq_len(max_k), function(k) DirichletMultinomial::dmn(count = counts, k = k, verbose = FALSE))
      criteria <- data.frame(k = seq_len(max_k), Laplace = vapply(fits, DirichletMultinomial::laplace, numeric(1)), AIC = vapply(fits, stats::AIC, numeric(1)), BIC = vapply(fits, stats::BIC, numeric(1)))
      best_k <- criteria$k[which.min(criteria$Laplace)]; best <- fits[[best_k]]
      assign <- apply(DirichletMultinomial::mixture(best), 1, which.max)
      assignments <- data.frame(Sample = rownames(counts), DMM_Cluster = factor(assign), stringsAsFactors = FALSE)
      list(ps_rank = x, counts = counts, criteria = criteria, best_k = best_k, best = best, assignments = assignments)
    }, ignoreInit = TRUE)
    fit_plot <- shiny::reactive({
      x <- fitted_models()$criteria; long <- tidyr::pivot_longer(x, -k, names_to = "Criterion", values_to = "Value")
      ggplot2::ggplot(long, ggplot2::aes(x = k, y = Value, colour = Criterion)) + ggplot2::geom_line(linewidth = 0.75) + ggplot2::geom_point(size = 2.5) +
        ggplot2::facet_wrap(~Criterion, scales = "free_y", ncol = 1) + ggplot2::scale_x_continuous(breaks = x$k) +
        ggplot2::labs(title = "DMM model selection", subtitle = paste("Best Laplace model: k =", fitted_models()$best_k), x = "Number of components", y = NULL) + zamp_theme() + ggplot2::theme(legend.position = "none")
    })
    driver_data <- shiny::reactive({
      fm <- fitted_models(); vals <- as.data.frame(DirichletMultinomial::fitted(fm$best)); vals$TaxonID <- rownames(vals)
      long <- tidyr::pivot_longer(vals, -TaxonID, names_to = "Component", values_to = "Weight"); long$Component <- sub("^.*?([0-9]+)$", "\\1", long$Component)
      tt <- as.data.frame(phyloseq::tax_table(fm$ps_rank), stringsAsFactors = FALSE); tt$TaxonID <- rownames(tt)
      if (identical(input$rank, "ASV") || !input$rank %in% names(tt)) tt$Taxon <- tt$TaxonID else {
        tt$Taxon <- as.character(tt[[input$rank]]); missing <- is.na(tt$Taxon) | !nzchar(tt$Taxon); tt$Taxon[missing] <- tt$TaxonID[missing]
      }
      long <- dplyr::left_join(long, tt[, c("TaxonID", "Taxon")], by = "TaxonID")
      long |> dplyr::group_by(Component) |> dplyr::slice_max(order_by = abs(Weight), n = 15, with_ties = FALSE) |> dplyr::ungroup()
    })
    drivers_plot <- shiny::reactive({
      df <- driver_data(); df$Taxon_display <- factor(paste(df$Taxon, df$Component, sep = "___"), levels = unique(paste(df$Taxon, df$Component, sep = "___")))
      ggplot2::ggplot(df, ggplot2::aes(x = Taxon_display, y = Weight, fill = factor(Component))) + ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
        ggplot2::coord_flip() + ggplot2::facet_wrap(~Component, scales = "free_y", ncol = 2) + ggplot2::scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
        ggplot2::scale_fill_manual(values = zamp_palette(length(unique(df$Component)))) + ggplot2::labs(title = "Top DMM community-type drivers", x = NULL, y = "Fitted component weight") + zamp_theme(10)
    })
    updated_ps <- shiny::reactive({
      fm <- fitted_models(); x <- ps(); meta <- as(phyloseq::sample_data(x), "data.frame"); meta$DMM_Cluster <- NA_character_
      meta[fm$assignments$Sample, "DMM_Cluster"] <- as.character(fm$assignments$DMM_Cluster); meta$DMM_Cluster <- factor(meta$DMM_Cluster); phyloseq::sample_data(x) <- phyloseq::sample_data(meta); x
    })
    output$fit_plot <- plotly::renderPlotly({ shiny::req(fit_plot()); zamp_plotly(fit_plot()) })
    output$drivers <- plotly::renderPlotly({ shiny::req(drivers_plot()); zamp_plotly(drivers_plot()) })
    output$best_model <- shiny::renderText({ shiny::req(fitted_models()); paste0("Selected model: k = ", fitted_models()$best_k, " (minimum Laplace criterion)") })
    output$assignments <- DT::renderDT({ shiny::req(fitted_models()); DT::datatable(fitted_models()$assignments, options = list(pageLength = 12), rownames = FALSE) })
    zamp_plot_download_server("fit_download", fit_plot, function() paste0("DMM_model_selection_", Sys.Date()))
    zamp_plot_download_server("drivers_download", drivers_plot, function() paste0("DMM_drivers_", Sys.Date()))
    output$assignments_download <- shiny::downloadHandler(filename = function() paste0("DMM_assignments_", Sys.Date(), ".tsv"), content = function(file) utils::write.table(fitted_models()$assignments, file, sep = "\t", row.names = FALSE, quote = FALSE))
    output$phyloseq_download <- shiny::downloadHandler(filename = function() paste0("phyloseq_with_DMM_clusters_", Sys.Date(), ".rds"), content = function(file) saveRDS(updated_ps(), file))
    shiny::observeEvent(input$run, { if (isTRUE(auto_save())) { shiny::req(fitted_models()); zamp_save_plot_to_dir(fit_plot(), output_dir(), "DMM_model_selection", "pdf", 6.5, 5); zamp_save_plot_to_dir(drivers_plot(), output_dir(), "DMM_drivers", "pdf", 10, 7); zamp_save_table_to_dir(fitted_models()$assignments, output_dir(), "DMM_assignments") } })
  })
}
