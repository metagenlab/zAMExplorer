mod_alpha_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Alpha-diversity settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("metrics"), "Metrics", choices = NULL, multiple = TRUE),
        shiny::selectInput(ns("group"), "Grouping variable", choices = NULL),
        shiny::actionButton(ns("generate"), "Calculate and plot", class = "btn-primary")
      )
    ),
    shiny::uiOutput(ns("richness_warning")),
    shiny::fluidRow(
      shinydashboard::box(title = "Alpha diversity", width = 12, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("plot"), height = "650px"),
        zamp_plot_download_ui(ns("download"), 8, 5.5)
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Pairwise Wilcoxon tests", width = 12, status = "primary", solidHeader = TRUE,
        DT::DTOutput(ns("stats")),
        shiny::downloadButton(ns("stats_download"), "Download statistics")
      )
    )
  )
}

zamp_alpha_calculate <- function(ps, measures) {
  measures <- unique(measures)
  warnings <- character()
  richness <- withCallingHandlers(
    phyloseq::estimate_richness(ps, measures = measures),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  richness <- as.data.frame(richness, check.names = FALSE)
  richness$Sample <- rownames(richness)

  available_metrics <- intersect(measures, colnames(richness))
  if (!length(available_metrics)) {
    stop("No requested alpha-diversity metrics could be calculated from this phyloseq object.")
  }

  meta <- zamp_metadata(ps)
  collisions <- intersect(available_metrics, colnames(meta))
  if (length(collisions)) {
    meta <- meta[, setdiff(colnames(meta), collisions), drop = FALSE]
  }

  table <- dplyr::left_join(richness, meta, by = "Sample")
  list(
    table = table,
    available_metrics = available_metrics,
    warnings = unique(warnings),
    replaced_metadata_metrics = collisions
  )
}

mod_alpha_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    alpha_result <- shiny::reactive({
      shiny::req(ps())
      tryCatch(
        zamp_alpha_calculate(
          ps(),
          measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
        ),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, paste("Alpha diversity could not be calculated:", conditionMessage(e))))
        }
      )
    })

    shiny::observeEvent(ps(), {
      shiny::req(ps())
      result <- alpha_result()
      metrics <- result$available_metrics
      groups <- zamp_categorical_variables(ps())
      selected_metrics <- intersect(c("Observed", "Shannon"), metrics)
      if (!length(selected_metrics) && length(metrics)) selected_metrics <- metrics[1]
      shiny::updateSelectInput(session, "metrics", choices = metrics, selected = selected_metrics)
      shiny::updateSelectInput(session, "group", choices = groups, selected = groups[1])
    }, ignoreInit = FALSE)

    output$richness_warning <- shiny::renderUI({
      result <- alpha_result()
      messages <- character()

      if (length(result$replaced_metadata_metrics)) {
        messages <- c(
          messages,
          paste0(
            "Existing metadata columns were replaced within the Alpha Diversity analysis by freshly calculated values: ",
            paste(result$replaced_metadata_metrics, collapse = ", "),
            ". The uploaded phyloseq object itself was not modified."
          )
        )
      }

      singleton_warning <- result$warnings[grepl("singleton|un-trimmed|richness estimates", result$warnings, ignore.case = TRUE)]
      if (length(singleton_warning)) {
        messages <- c(
          messages,
          paste0(
            "Richness warning from phyloseq: this dataset contains no singletons. Chao1/ACE can be unreliable when low-abundance taxa were removed before analysis. Observed, Shannon, Simpson and InvSimpson are still available."
          )
        )
      }

      if (!length(messages)) return(NULL)
      shiny::fluidRow(
        shinydashboard::box(
          title = "Alpha-diversity note",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          lapply(unique(messages), shiny::tags$p)
        )
      )
    })

    plot_data <- shiny::eventReactive(input$generate, {
      shiny::req(input$metrics, input$group)
      result <- alpha_result()
      df <- result$table
      missing_metrics <- setdiff(input$metrics, colnames(df))
      shiny::validate(shiny::need(!length(missing_metrics), paste("These selected alpha metrics could not be calculated:", paste(missing_metrics, collapse = ", "))))
      shiny::validate(shiny::need(input$group %in% colnames(df), "The selected grouping variable is not available in the current phyloseq metadata."))

      long <- tidyr::pivot_longer(
        df,
        cols = dplyr::all_of(input$metrics),
        names_to = "Metric",
        values_to = "Diversity"
      )
      long$Group <- factor(long[[input$group]])
      long <- long[!is.na(long$Group) & is.finite(long$Diversity), , drop = FALSE]
      long$Group <- droplevels(long$Group)
      shiny::validate(shiny::need(nrow(long) > 0, "No finite alpha-diversity values remain for the selected settings."))
      long
    }, ignoreInit = TRUE)

    alpha_plot <- shiny::reactive({
      df <- plot_data()
      pal <- stats::setNames(zamp_palette(nlevels(df$Group)), levels(df$Group))
      ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Diversity, colour = Group)) +
        ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.55) +
        ggplot2::geom_jitter(width = 0.14, size = 2.1, alpha = 0.72) +
        ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
        ggplot2::facet_wrap(~Metric, scales = "free_y", nrow = 1) +
        ggplot2::labs(title = "Alpha diversity", subtitle = paste("Grouped by", input$group), x = NULL, y = NULL) +
        zamp_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
    })

    stats_df <- shiny::reactive({
      df <- plot_data()
      out <- lapply(split(df, df$Metric), function(d) {
        d$Group <- droplevels(d$Group)
        groups <- levels(d$Group)
        if (length(groups) < 2) return(NULL)
        comps <- utils::combn(groups, 2, simplify = FALSE)
        tab <- lapply(comps, function(cp) {
          x <- d$Diversity[d$Group == cp[1]]
          y <- d$Diversity[d$Group == cp[2]]
          wt <- tryCatch(stats::wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
          data.frame(
            Group_1 = cp[1],
            Group_2 = cp[2],
            n_1 = length(x),
            n_2 = length(y),
            p_value = if (is.null(wt)) NA_real_ else wt$p.value,
            stringsAsFactors = FALSE
          )
        })
        ans <- dplyr::bind_rows(tab)
        if (!nrow(ans)) return(NULL)
        ans$Metric <- unique(as.character(d$Metric))[1]
        ans$p_adjusted_BH <- stats::p.adjust(ans$p_value, method = "BH")
        ans[, c("Metric", "Group_1", "Group_2", "n_1", "n_2", "p_value", "p_adjusted_BH")]
      })
      ans <- dplyr::bind_rows(out)
      if (!nrow(ans)) {
        ans <- data.frame(
          Metric = character(), Group_1 = character(), Group_2 = character(),
          n_1 = integer(), n_2 = integer(), p_value = numeric(), p_adjusted_BH = numeric(),
          stringsAsFactors = FALSE
        )
      }
      ans
    })

    output$plot <- plotly::renderPlotly({
      shiny::req(alpha_plot())
      tryCatch(
        zamp_plotly(alpha_plot()),
        error = function(e) shiny::validate(shiny::need(FALSE, paste("Alpha-diversity plot could not be rendered:", conditionMessage(e))))
      )
    })

    output$stats <- DT::renderDT({
      x <- stats_df()
      if (!nrow(x)) {
        return(DT::datatable(
          data.frame(Message = "Pairwise Wilcoxon tests require at least two non-missing groups."),
          options = list(dom = "t"), rownames = FALSE
        ))
      }
      x$p_value <- signif(x$p_value, 4)
      x$p_adjusted_BH <- signif(x$p_adjusted_BH, 4)
      DT::datatable(x, options = list(pageLength = 15), rownames = FALSE)
    })

    zamp_plot_download_server("download", alpha_plot, function() paste0("alpha_diversity_", Sys.Date()))

    output$stats_download <- shiny::downloadHandler(
      filename = function() paste0("alpha_diversity_stats_", Sys.Date(), ".tsv"),
      content = function(file) {
        x <- stats_df()
        utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    shiny::observeEvent(input$generate, {
      if (isTRUE(auto_save())) {
        tryCatch({
          shiny::req(alpha_plot())
          pfile <- zamp_save_plot_to_dir(alpha_plot(), output_dir(), "alpha_diversity", "pdf", 8, 5.5)
          tfile <- zamp_save_table_to_dir(stats_df(), output_dir(), "alpha_diversity_stats")
          shiny::showNotification(paste("Saved:", basename(pfile), "and", basename(tfile)), type = "message")
        }, error = function(e) {
          shiny::showNotification(paste("Alpha Diversity ran, but automatic saving failed:", conditionMessage(e)), type = "warning", duration = NULL)
        })
      }
    }, ignoreInit = TRUE)
  })
}
