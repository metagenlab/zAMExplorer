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

mod_alpha_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    alpha_table <- shiny::reactive({
      shiny::req(ps())
      richness <- phyloseq::estimate_richness(ps(), measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
      richness$Sample <- rownames(richness)
      dplyr::left_join(richness, zamp_metadata(ps()), by = "Sample")
    })

    shiny::observeEvent(ps(), {
      shiny::req(ps())
      metrics <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
      groups <- zamp_categorical_variables(ps())
      shiny::updateSelectInput(session, "metrics", choices = metrics, selected = intersect(c("Observed", "Shannon"), metrics))
      shiny::updateSelectInput(session, "group", choices = groups, selected = groups[1])
    })

    plot_data <- shiny::eventReactive(input$generate, {
      shiny::req(input$metrics, input$group)
      df <- alpha_table()
      long <- tidyr::pivot_longer(df, dplyr::all_of(input$metrics), names_to = "Metric", values_to = "Diversity")
      long$Group <- factor(long[[input$group]])
      long <- long[!is.na(long$Group) & is.finite(long$Diversity), , drop = FALSE]
      shiny::validate(shiny::need(nlevels(long$Group) >= 2, "The selected grouping variable needs at least two non-missing groups."))
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
        groups <- levels(droplevels(d$Group))
        if (length(groups) < 2) return(NULL)
        comps <- utils::combn(groups, 2, simplify = FALSE)
        tab <- lapply(comps, function(cp) {
          x <- d$Diversity[d$Group == cp[1]]
          y <- d$Diversity[d$Group == cp[2]]
          wt <- tryCatch(stats::wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
          data.frame(Group_1 = cp[1], Group_2 = cp[2], n_1 = length(x), n_2 = length(y), p_value = if (is.null(wt)) NA_real_ else wt$p.value)
        })
        ans <- dplyr::bind_rows(tab)
        ans$Metric <- unique(d$Metric)
        ans$p_adjusted_BH <- stats::p.adjust(ans$p_value, method = "BH")
        ans[, c("Metric", "Group_1", "Group_2", "n_1", "n_2", "p_value", "p_adjusted_BH")]
      })
      dplyr::bind_rows(out)
    })

    output$plot <- plotly::renderPlotly({ shiny::req(alpha_plot()); zamp_plotly(alpha_plot()) })
    output$stats <- DT::renderDT({
      x <- stats_df(); x$p_value <- signif(x$p_value, 4); x$p_adjusted_BH <- signif(x$p_adjusted_BH, 4)
      DT::datatable(x, options = list(pageLength = 15), rownames = FALSE)
    })

    zamp_plot_download_server("download", alpha_plot, function() paste0("alpha_diversity_", Sys.Date()))
    output$stats_download <- shiny::downloadHandler(
      filename = function() paste0("alpha_diversity_stats_", Sys.Date(), ".tsv"),
      content = function(file) utils::write.table(stats_df(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    )

    shiny::observeEvent(input$generate, {
      if (isTRUE(auto_save())) {
        shiny::req(alpha_plot())
        pfile <- zamp_save_plot_to_dir(alpha_plot(), output_dir(), "alpha_diversity", "pdf", 8, 5.5)
        tfile <- zamp_save_table_to_dir(stats_df(), output_dir(), "alpha_diversity_stats")
        shiny::showNotification(paste("Saved:", basename(pfile), "and", basename(tfile)), type = "message")
      }
    })
  })
}
