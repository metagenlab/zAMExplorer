mod_maaslin2_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "MaAsLin2 settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::selectInput(ns("fixed"), "Fixed effects", choices = NULL, multiple = TRUE),
        shiny::selectInput(ns("random"), "Random effects (optional)", choices = NULL, multiple = TRUE),
        shiny::textInput(ns("reference"), "Reference level(s), optional", placeholder = "group,Control;sex,Female"),
        shiny::selectInput(ns("normalization"), "Normalization", c("TSS", "CLR", "CSS", "NONE", "TMM"), selected = "TSS"),
        shiny::selectInput(ns("transform"), "Transform", c("LOG", "LOGIT", "AST", "NONE"), selected = "LOG"),
        shiny::selectInput(ns("method"), "Analysis method", c("Linear model" = "LM", "Compound Poisson" = "CPLM", "Negative binomial" = "NEGBIN", "Zero-inflated negative binomial" = "ZINB"), selected = "LM"),
        shiny::numericInput(ns("min_abundance"), "Minimum abundance", 0, min = 0, step = 0.0001),
        shiny::numericInput(ns("min_prevalence"), "Minimum prevalence (%)", 10, min = 0, max = 100, step = 1),
        shiny::numericInput(ns("q_threshold"), "q-value threshold", 0.10, min = 0.001, max = 1, step = 0.01),
        shiny::checkboxInput(ns("standardize"), "Standardize continuous metadata", TRUE),
        shiny::actionButton(ns("run"), "Run MaAsLin2", class = "btn-primary")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Association coefficients", width = 7, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("coef_plot"), height = "700px"), zamp_plot_download_ui(ns("coef_download"), 8, 7)),
      shinydashboard::box(title = "MaAsLin2 results", width = 5, status = "primary", solidHeader = TRUE,
        shiny::verbatimTextOutput(ns("summary")), DT::DTOutput(ns("results")), shiny::downloadButton(ns("results_download"), "Download all results"))
    )
  )
}

zamp_validate_reference <- function(reference, metadata) {
  reference <- trimws(reference %||% "")
  if (!nzchar(reference)) return(NULL)
  entries <- trimws(strsplit(reference, ";", fixed = TRUE)[[1]])
  entries <- entries[nzchar(entries)]
  bad <- character()
  for (entry in entries) {
    parts <- trimws(strsplit(entry, ",", fixed = TRUE)[[1]])
    if (length(parts) != 2 || !nzchar(parts[1]) || !nzchar(parts[2])) {
      bad <- c(bad, entry)
      next
    }
    variable <- parts[1]
    level <- parts[2]
    if (!variable %in% colnames(metadata) || !level %in% as.character(unique(metadata[[variable]]))) {
      bad <- c(bad, entry)
    }
  }
  if (length(bad)) {
    stop("Invalid reference specification: ", paste(bad, collapse = "; "), ". Use variable,referenceLevel;variable,referenceLevel.")
  }
  paste(entries, collapse = ";")
}

mod_maaslin2_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps()))
      vars <- zamp_all_metadata_variables(ps())
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = if ("Genus" %in% ranks) "Genus" else ranks[1])
      shiny::updateSelectInput(session, "fixed", choices = vars)
      shiny::updateSelectInput(session, "random", choices = vars)
    })

    fit <- shiny::eventReactive(input$run, {
      shiny::req(ps(), input$rank, input$fixed)
      tryCatch({
        shiny::validate(shiny::need(length(input$fixed) >= 1, "Select at least one fixed effect."))
        shiny::validate(shiny::need(length(intersect(input$fixed, input$random)) == 0, "A variable cannot be both a fixed and random effect."))

        x <- zamp_aggregate_rank(ps(), input$rank)
        data <- zamp_otu_matrix(x, samples_in_rows = TRUE)
        labels <- zamp_taxon_labels(x, input$rank)
        label_map <- stats::setNames(labels$Taxon, labels$TaxonID)
        colnames(data) <- make.unique(unname(label_map[colnames(data)]), sep = "_")

        metadata <- zamp_metadata(x)
        rownames(metadata) <- metadata$Sample
        metadata <- metadata[rownames(data), , drop = FALSE]
        metadata$Sample <- NULL

        selected_effects <- unique(c(input$fixed, input$random))
        missing_effects <- setdiff(selected_effects, colnames(metadata))
        shiny::validate(shiny::need(!length(missing_effects), paste("Selected metadata variables are missing:", paste(missing_effects, collapse = ", "))))

        keep <- stats::complete.cases(metadata[, selected_effects, drop = FALSE])
        data <- data[keep, , drop = FALSE]
        metadata <- metadata[keep, , drop = FALSE]
        shiny::validate(shiny::need(nrow(data) >= 4, "At least four samples with complete selected metadata are required for MaAsLin2."))
        shiny::validate(shiny::need(ncol(data) >= 1, "No taxa remain for MaAsLin2."))

        for (v in input$fixed) {
          values <- metadata[[v]]
          if (is.character(values) || is.factor(values) || is.logical(values)) {
            shiny::validate(shiny::need(length(unique(values)) >= 2, paste("Fixed effect", v, "has fewer than two levels after removing missing samples.")))
          } else {
            shiny::validate(shiny::need(stats::sd(values, na.rm = TRUE) > 0, paste("Fixed effect", v, "has no variation after removing missing samples.")))
          }
        }
        if (length(input$random)) {
          for (v in input$random) {
            shiny::validate(shiny::need(length(unique(metadata[[v]])) >= 2, paste("Random effect", v, "has fewer than two levels after removing missing samples.")))
          }
        }

        ref <- zamp_validate_reference(input$reference, metadata)
        random <- if (length(input$random)) input$random else NULL
        run_dir <- if (isTRUE(auto_save())) {
          file.path(output_dir(), paste0("MaAsLin2_", base::format(Sys.time(), "%Y%m%d_%H%M%S")))
        } else {
          tempfile("zAMPExplorer_MaAsLin2_")
        }
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

        fit_obj <- Maaslin2::Maaslin2(
          input_data = data,
          input_metadata = metadata,
          output = run_dir,
          min_abundance = input$min_abundance,
          min_prevalence = input$min_prevalence / 100,
          normalization = input$normalization,
          transform = input$transform,
          analysis_method = input$method,
          max_significance = input$q_threshold,
          random_effects = random,
          fixed_effects = input$fixed,
          correction = "BH",
          standardize = isTRUE(input$standardize),
          cores = 1,
          plot_heatmap = FALSE,
          plot_scatter = FALSE,
          save_scatter = FALSE,
          save_models = FALSE,
          reference = ref
        )

        results <- fit_obj$results
        if (is.null(results)) {
          all_results <- file.path(run_dir, "all_results.tsv")
          shiny::validate(shiny::need(file.exists(all_results), "MaAsLin2 did not return a results table."))
          results <- utils::read.delim(all_results, check.names = FALSE)
        }
        results <- as.data.frame(results, stringsAsFactors = FALSE)
        needed_columns <- c("feature", "metadata", "value", "coef", "stderr", "pval", "qval")
        missing_columns <- setdiff(needed_columns, colnames(results))
        shiny::validate(shiny::need(!length(missing_columns), paste("MaAsLin2 results are missing expected columns:", paste(missing_columns, collapse = ", "))))
        results <- results[order(results$qval, results$pval), , drop = FALSE]
        list(results = results, run_dir = run_dir, n_used = nrow(data), n_removed = sum(!keep))
      }, error = function(e) {
        shiny::validate(shiny::need(FALSE, paste("MaAsLin2 failed:", conditionMessage(e))))
      })
    }, ignoreInit = TRUE)

    coef_data <- shiny::reactive({
      results <- fit()$results
      shiny::validate(shiny::need(nrow(results) > 0, "No MaAsLin2 associations were returned."))
      n <- min(40, nrow(results))
      x <- results[seq_len(n), , drop = FALSE]
      x$Association <- ifelse(is.na(x$value) | !nzchar(as.character(x$value)), x$metadata, paste0(x$metadata, ": ", x$value))
      x$FeatureAssociation <- paste(x$feature, x$Association, sep = "___")
      x$Significant <- x$qval <= input$q_threshold
      x
    })

    coef_plot <- shiny::reactive({
      x <- coef_data()
      x$FeatureAssociation <- factor(x$FeatureAssociation, levels = rev(x$FeatureAssociation))
      ggplot2::ggplot(
        x,
        ggplot2::aes(
          x = coef, y = FeatureAssociation, colour = Significant,
          text = paste0("Feature: ", feature, "<br>Metadata: ", Association, "<br>Coefficient: ", signif(coef, 4), "<br>q-value: ", signif(qval, 3))
        )
      ) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey45") +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = coef - 1.96 * stderr, xmax = coef + 1.96 * stderr), orientation = "y", width = 0.18, linewidth = 0.45) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_y_discrete(labels = function(v) sub("___", "  |  ", v, fixed = TRUE)) +
        ggplot2::scale_colour_manual(values = c("FALSE" = "grey55", "TRUE" = "#0072B2"), labels = c("q above threshold", "q at/below threshold")) +
        ggplot2::labs(title = "MaAsLin2 association coefficients", subtitle = paste0("Top ", nrow(x), " associations ordered by q-value; 95% model CI"), x = "Model coefficient", y = NULL, colour = NULL) +
        zamp_theme(10) +
        ggplot2::theme(legend.position = "top")
    })

    output$coef_plot <- plotly::renderPlotly({
      shiny::req(coef_plot())
      tryCatch(
        zamp_plotly(coef_plot(), tooltip = "text"),
        error = function(e) shiny::validate(shiny::need(FALSE, paste("MaAsLin2 coefficient plot could not be rendered:", conditionMessage(e))))
      )
    })
    output$results <- DT::renderDT({
      shiny::req(fit())
      DT::datatable(fit()$results, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE, filter = "top")
    })
    output$summary <- shiny::renderText({
      r <- fit()$results
      sprintf(
        "%d tested associations\n%d associations with q <= %.3g\nSamples used: %d; excluded for missing selected metadata: %d\nMaAsLin2 output: %s",
        nrow(r), sum(r$qval <= input$q_threshold, na.rm = TRUE), input$q_threshold,
        fit()$n_used, fit()$n_removed, fit()$run_dir
      )
    })

    zamp_plot_download_server("coef_download", coef_plot, function() paste0("MaAsLin2_coefficients_", Sys.Date()))
    output$results_download <- shiny::downloadHandler(
      filename = function() paste0("MaAsLin2_all_results_", Sys.Date(), ".tsv"),
      content = function(file) {
        shiny::req(fit())
        utils::write.table(fit()$results, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    shiny::observeEvent(input$run, {
      if (isTRUE(auto_save())) {
        tryCatch({
          shiny::req(fit())
          zamp_save_plot_to_dir(coef_plot(), fit()$run_dir, "zAMPExplorer_MaAsLin2_coefficients", "pdf", 8, 7)
          zamp_save_table_to_dir(fit()$results, fit()$run_dir, "zAMPExplorer_MaAsLin2_results")
        }, error = function(e) {
          shiny::showNotification(paste("MaAsLin2 ran, but zAMPExplorer automatic export failed:", conditionMessage(e)), type = "warning", duration = NULL)
        })
      }
    }, ignoreInit = TRUE)
  })
}
