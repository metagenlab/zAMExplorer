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
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps()))
      vars <- zamp_all_metadata_variables(ps())
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = if ("Genus" %in% ranks) "Genus" else ranks[1])
      shiny::updateSelectInput(session, "explanatory", choices = vars, selected = vars[seq_len(min(2, length(vars)))])
      shiny::updateSelectInput(session, "colour", choices = vars, selected = vars[1])
    })

    result <- shiny::eventReactive(input$run, {
      shiny::req(ps(), input$rank, input$explanatory, input$colour)
      tryCatch({
        x <- zamp_aggregate_rank(ps(), input$rank)
        comm <- zamp_otu_matrix(x, samples_in_rows = TRUE)
        meta <- zamp_metadata(x)
        rownames(meta) <- meta$Sample
        meta <- meta[rownames(comm), , drop = FALSE]

        missing_vars <- setdiff(unique(c(input$explanatory, input$colour)), colnames(meta))
        shiny::validate(shiny::need(!length(missing_vars), paste("Selected metadata variables are missing:", paste(missing_vars, collapse = ", "))))

        needed <- unique(c(input$explanatory, input$colour))
        keep <- stats::complete.cases(meta[, needed, drop = FALSE])
        comm <- comm[keep, , drop = FALSE]
        meta <- meta[keep, , drop = FALSE]
        shiny::validate(shiny::need(nrow(comm) >= 4, "Too few complete samples remain for constrained ordination."))

        env <- meta[, input$explanatory, drop = FALSE]
        bad_variables <- vapply(env, function(v) {
          if (is.numeric(v)) {
            s <- stats::sd(v, na.rm = TRUE)
            !is.finite(s) || s == 0
          } else {
            length(unique(v)) < 2
          }
        }, logical(1))
        shiny::validate(shiny::need(!any(bad_variables), paste("These explanatory variables have no variation after excluding missing samples:", paste(names(env)[bad_variables], collapse = ", "))))

        names(env) <- make.names(names(env), unique = TRUE)
        design <- stats::model.matrix(~ ., data = env)
        design_rank <- qr(design)$rank
        predictor_df <- max(0L, design_rank - 1L)
        shiny::validate(shiny::need(predictor_df >= 1, "The selected explanatory variables do not produce an estimable constrained model."))
        shiny::validate(shiny::need(predictor_df < nrow(comm) - 1L, paste0("The selected explanatory variables use ", predictor_df, " model degrees of freedom but only ", nrow(comm), " complete samples remain. Select fewer/simpler explanatory variables.")))

        comm <- zamp_transform_matrix(comm, input$transform)
        shiny::validate(
          shiny::need(ncol(comm) >= 2, "At least two taxa/features are required for constrained ordination."),
          shiny::need(all(is.finite(comm)), "The transformed community matrix contains non-finite values.")
        )

        formula <- stats::as.formula("comm ~ .")
        fit <- if (identical(input$method, "RDA")) {
          vegan::rda(formula, data = env)
        } else {
          vegan::capscale(formula, data = env, distance = input$distance, add = TRUE)
        }

        sites <- as.data.frame(vegan::scores(fit, display = "sites", choices = 1:2, scaling = 2))
        shiny::validate(shiny::need(ncol(sites) >= 2, "The constrained ordination returned fewer than two site axes."))
        sites$Sample <- rownames(sites)
        axis_names <- names(sites)[1:2]
        sites <- dplyr::left_join(sites, meta, by = "Sample")
        sites$Colour <- factor(sites[[input$colour]])

        eig <- tryCatch(vegan::eigenvals(fit, model = "constrained"), error = function(e) numeric())
        pct <- if (length(eig) >= 2 && sum(eig) > 0) eig[1:2] / sum(eig) * 100 else c(NA_real_, NA_real_)
        xlab <- if (is.finite(pct[1])) sprintf("Axis 1 (%.1f%% constrained)", pct[1]) else axis_names[1]
        ylab <- if (is.finite(pct[2])) sprintf("Axis 2 (%.1f%% constrained)", pct[2]) else axis_names[2]

        species <- tryCatch(as.data.frame(vegan::scores(fit, display = "species", choices = 1:2, scaling = 2)), error = function(e) NULL)
        if (!is.null(species) && nrow(species) && input$top_taxa > 0) {
          species$TaxonID <- rownames(species)
          score_cols <- names(species)[1:2]
          species$length <- sqrt(species[[score_cols[1]]]^2 + species[[score_cols[2]]]^2)
          species <- utils::head(species[order(species$length, decreasing = TRUE), , drop = FALSE], input$top_taxa)
          tt <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = FALSE)
          tt$TaxonID <- rownames(tt)
          tt$Taxon <- if (identical(input$rank, "ASV") || !input$rank %in% names(tt)) tt$TaxonID else as.character(tt[[input$rank]])
          missing <- is.na(tt$Taxon) | !nzchar(tt$Taxon)
          tt$Taxon[missing] <- tt$TaxonID[missing]
          species <- dplyr::left_join(species, tt[, c("TaxonID", "Taxon")], by = "TaxonID")
        } else {
          species <- NULL
        }

        pal <- stats::setNames(zamp_palette(nlevels(sites$Colour)), levels(sites$Colour))
        p <- ggplot2::ggplot(sites, ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]], colour = Colour, text = Sample)) +
          ggplot2::geom_point(size = 3.1, alpha = 0.9) +
          ggplot2::scale_colour_manual(values = pal) +
          ggplot2::labs(
            title = if (identical(input$method, "RDA")) "Redundancy analysis (RDA)" else "Distance-based RDA (dbRDA)",
            subtitle = paste("Constraints:", paste(input$explanatory, collapse = " + ")),
            x = xlab, y = ylab, colour = input$colour
          ) +
          zamp_theme() +
          ggplot2::theme(aspect.ratio = 1)

        if (!is.null(species) && nrow(species)) {
          score_cols <- names(species)[1:2]
          p <- p +
            ggplot2::geom_segment(
              data = species,
              ggplot2::aes(x = 0, y = 0, xend = .data[[score_cols[1]]], yend = .data[[score_cols[2]]]),
              inherit.aes = FALSE,
              arrow = grid::arrow(length = grid::unit(0.16, "cm")),
              linewidth = 0.45,
              colour = "grey30"
            ) +
            ggrepel::geom_text_repel(
              data = species,
              ggplot2::aes(x = .data[[score_cols[1]]], y = .data[[score_cols[2]]], label = Taxon),
              inherit.aes = FALSE,
              size = 3, max.overlaps = 30, colour = "grey15"
            )
        }

        overall <- stats::anova(fit, permutations = input$permutations)
        terms <- stats::anova(fit, permutations = input$permutations, by = "margin")
        anova_df <- dplyr::bind_rows(
          data.frame(Test = "Overall", Term = rownames(overall), as.data.frame(overall), row.names = NULL, check.names = FALSE),
          data.frame(Test = "Marginal", Term = rownames(terms), as.data.frame(terms), row.names = NULL, check.names = FALSE)
        )
        list(plot = p, anova = anova_df, n_used = nrow(comm), n_removed = sum(!keep), predictor_df = predictor_df)
      }, error = function(e) {
        shiny::validate(shiny::need(FALSE, paste("RDA/dbRDA analysis failed:", conditionMessage(e))))
      })
    }, ignoreInit = TRUE)

    rda_plot <- shiny::reactive({
      shiny::req(result())
      result()$plot
    })
    output$plot <- plotly::renderPlotly({
      shiny::req(rda_plot())
      tryCatch(
        zamp_plotly(rda_plot(), tooltip = c("text", "x", "y", "colour")),
        error = function(e) shiny::validate(shiny::need(FALSE, paste("RDA/dbRDA plot could not be rendered:", conditionMessage(e))))
      )
    })
    output$note <- shiny::renderText({
      shiny::req(result())
      sprintf(
        "Samples used: %d; samples excluded for missing selected metadata: %d; constrained model degrees of freedom: %d",
        result()$n_used, result()$n_removed, result()$predictor_df
      )
    })
    output$anova <- DT::renderDT({
      shiny::req(result())
      DT::datatable(result()$anova, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
    })

    zamp_plot_download_server("download", rda_plot, function() paste0(input$method, "_ordination_", Sys.Date()))
    output$anova_download <- shiny::downloadHandler(
      filename = function() paste0(input$method, "_permutation_tests_", Sys.Date(), ".tsv"),
      content = function(file) {
        shiny::req(result())
        utils::write.table(result()$anova, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    shiny::observeEvent(input$run, {
      if (isTRUE(auto_save())) {
        tryCatch({
          shiny::req(result())
          zamp_save_plot_to_dir(result()$plot, output_dir(), paste0(input$method, "_ordination"), "pdf", 8, 6.5)
          zamp_save_table_to_dir(result()$anova, output_dir(), paste0(input$method, "_permutation_tests"))
        }, error = function(e) {
          shiny::showNotification(paste("RDA/dbRDA ran, but automatic export failed:", conditionMessage(e)), type = "warning", duration = NULL)
        })
      }
    }, ignoreInit = TRUE)
  })
}
