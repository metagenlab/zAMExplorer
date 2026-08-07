mod_rda_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Constrained ordination settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("method"), "Method", c("RDA" = "RDA", "distance-based RDA (dbRDA)" = "dbRDA")),
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::selectInput(ns("transform"), "Community transformation", c("Hellinger" = "hellinger", "Compositional" = "compositional", "CLR" = "clr", "Identity" = "identity")),
        shiny::selectInput(ns("distance"), "dbRDA distance", c("Bray-Curtis" = "bray", "Jaccard" = "jaccard", "Euclidean" = "euclidean")),
        shiny::selectInput(ns("explanatory"), "Explanatory metadata variables", choices = NULL, multiple = TRUE),
        shiny::selectInput(ns("colour"), "Colour samples by", choices = NULL),
        shiny::numericInput(ns("top_taxa"), "Number of taxa arrows", 12, min = 0, max = 50),
        shiny::numericInput(ns("permutations"), "Permutations", 999, min = 99, max = 99999, step = 100),
        shiny::actionButton(ns("run"), "Run constrained ordination", class = "btn-primary")
      )
    ),
    shiny::fluidRow(shinydashboard::box(title = "RDA / dbRDA", width = 12, status = "primary", solidHeader = TRUE,
      plotly::plotlyOutput(ns("plot"), height = "700px"), shiny::verbatimTextOutput(ns("note")), zamp_plot_download_ui(ns("download"), 8, 6.5))),
    shiny::fluidRow(shinydashboard::box(title = "Permutation tests", width = 12, status = "primary", solidHeader = TRUE,
      DT::DTOutput(ns("anova")), shiny::downloadButton(ns("anova_download"), "Download permutation tests")))
  )
}

mod_rda_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps()); ranks <- c("ASV", phyloseq::rank_names(ps())); vars <- zamp_all_metadata_variables(ps())
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = if ("Genus" %in% ranks) "Genus" else ranks[1])
      shiny::updateSelectInput(session, "explanatory", choices = vars, selected = vars[seq_len(min(2, length(vars)))])
      shiny::updateSelectInput(session, "colour", choices = vars, selected = vars[1])
    })
    result <- shiny::eventReactive(input$run, {
      shiny::req(ps(), input$rank, input$explanatory, input$colour)
      x <- zamp_aggregate_rank(ps(), input$rank); comm <- zamp_otu_matrix(x, samples_in_rows = TRUE); comm <- zamp_transform_matrix(comm, input$transform)
      meta <- zamp_metadata(x); rownames(meta) <- meta$Sample; meta <- meta[rownames(comm), , drop = FALSE]
      needed <- unique(c(input$explanatory, input$colour)); keep <- stats::complete.cases(meta[, needed, drop = FALSE]); comm <- comm[keep, , drop = FALSE]; meta <- meta[keep, , drop = FALSE]
      shiny::validate(shiny::need(nrow(comm) >= 4, "Too few complete samples remain for constrained ordination."))
      shiny::validate(shiny::need(length(input$explanatory) < nrow(comm) - 1, "Too many explanatory variables for the available samples."))
      env <- meta[, input$explanatory, drop = FALSE]; names(env) <- make.names(names(env), unique = TRUE); formula <- stats::as.formula("comm ~ .")
      if (identical(input$method, "RDA")) fit <- vegan::rda(formula, data = env) else fit <- vegan::capscale(formula, data = env, distance = input$distance, add = TRUE)
      sites <- as.data.frame(vegan::scores(fit, display = "sites", choices = 1:2, scaling = 2)); sites$Sample <- rownames(sites); axis_names <- names(sites)[1:2]
      sites <- dplyr::left_join(sites, meta, by = "Sample"); sites$Colour <- factor(sites[[input$colour]])
      eig <- tryCatch(vegan::eigenvals(fit, model = "constrained"), error = function(e) numeric()); pct <- if (length(eig) >= 2 && sum(eig) > 0) eig[1:2] / sum(eig) * 100 else c(NA_real_, NA_real_)
      xlab <- if (is.finite(pct[1])) sprintf("Axis 1 (%.1f%% constrained)", pct[1]) else axis_names[1]; ylab <- if (is.finite(pct[2])) sprintf("Axis 2 (%.1f%% constrained)", pct[2]) else axis_names[2]
      species <- tryCatch(as.data.frame(vegan::scores(fit, display = "species", choices = 1:2, scaling = 2)), error = function(e) NULL)
      if (!is.null(species) && nrow(species)) {
        species$TaxonID <- rownames(species); score_cols <- names(species)[1:2]; species$length <- sqrt(species[[score_cols[1]]]^2 + species[[score_cols[2]]]^2)
        species <- utils::head(species[order(species$length, decreasing = TRUE), , drop = FALSE], input$top_taxa)
        tt <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = FALSE); tt$TaxonID <- rownames(tt); tt$Taxon <- if (identical(input$rank, "ASV") || !input$rank %in% names(tt)) tt$TaxonID else as.character(tt[[input$rank]])
        missing <- is.na(tt$Taxon) | !nzchar(tt$Taxon); tt$Taxon[missing] <- tt$TaxonID[missing]; species <- dplyr::left_join(species, tt[, c("TaxonID", "Taxon")], by = "TaxonID")
      }
      pal <- stats::setNames(zamp_palette(nlevels(sites$Colour)), levels(sites$Colour))
      p <- ggplot2::ggplot(sites, ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]], colour = Colour, text = Sample)) +
        ggplot2::geom_point(size = 3.1, alpha = 0.9) + ggplot2::scale_colour_manual(values = pal) +
        ggplot2::labs(title = if (identical(input$method, "RDA")) "Redundancy analysis (RDA)" else "Distance-based RDA (dbRDA)", subtitle = paste("Constraints:", paste(input$explanatory, collapse = " + ")), x = xlab, y = ylab, colour = input$colour) + zamp_theme() + ggplot2::theme(aspect.ratio = 1)
      if (!is.null(species) && nrow(species)) {
        score_cols <- names(species)[1:2]
        p <- p + ggplot2::geom_segment(data = species, ggplot2::aes(x = 0, y = 0, xend = .data[[score_cols[1]]], yend = .data[[score_cols[2]]]), inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.16, "cm")), linewidth = 0.45, colour = "grey30") +
          ggrepel::geom_text_repel(data = species, ggplot2::aes(x = .data[[score_cols[1]]], y = .data[[score_cols[2]]], label = Taxon), inherit.aes = FALSE, size = 3, max.overlaps = 30, colour = "grey15")
      }
      overall <- stats::anova(fit, permutations = input$permutations); terms <- stats::anova(fit, permutations = input$permutations, by = "margin")
      anova_df <- dplyr::bind_rows(data.frame(Test = "Overall", Term = rownames(overall), as.data.frame(overall), row.names = NULL, check.names = FALSE), data.frame(Test = "Marginal", Term = rownames(terms), as.data.frame(terms), row.names = NULL, check.names = FALSE))
      list(plot = p, anova = anova_df, n_used = nrow(comm), n_removed = sum(!keep))
    }, ignoreInit = TRUE)
    rda_plot <- shiny::reactive({ shiny::req(result()); result()$plot })
    output$plot <- plotly::renderPlotly({ shiny::req(rda_plot()); zamp_plotly(rda_plot(), tooltip = c("text", "x", "y", "colour")) })
    output$note <- shiny::renderText({ shiny::req(result()); sprintf("Samples used: %d; samples excluded for missing selected metadata: %d", result()$n_used, result()$n_removed) })
    output$anova <- DT::renderDT({ shiny::req(result()); DT::datatable(result()$anova, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE) })
    zamp_plot_download_server("download", rda_plot, function() paste0(input$method, "_ordination_", Sys.Date()))
    output$anova_download <- shiny::downloadHandler(filename = function() paste0(input$method, "_permutation_tests_", Sys.Date(), ".tsv"), content = function(file) utils::write.table(result()$anova, file, sep = "\t", row.names = FALSE, quote = FALSE))
    shiny::observeEvent(input$run, { if (isTRUE(auto_save())) { shiny::req(result()); zamp_save_plot_to_dir(result()$plot, output_dir(), paste0(input$method, "_ordination"), "pdf", 8, 6.5); zamp_save_table_to_dir(result()$anova, output_dir(), paste0(input$method, "_permutation_tests")) } })
  })
}
