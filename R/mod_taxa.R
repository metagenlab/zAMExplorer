mod_taxa_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Taxa overview settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::selectInput(ns("group"), "Grouping variable", choices = NULL),
        shiny::numericInput(ns("detection"), "Detection threshold (%)", 0.1, min = 0, max = 100, step = 0.05),
        shiny::numericInput(ns("prevalence"), "Core prevalence threshold (%)", 20, min = 0, max = 100, step = 1),
        shiny::actionButton(ns("generate"), "Calculate taxa overview", class = "btn-primary")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Prevalence and abundance", width = 7, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("plot"), height = "560px"),
        zamp_plot_download_ui(ns("download"), 7, 5.5)
      ),
      shinydashboard::box(title = "Taxa prevalence table", width = 5, status = "primary", solidHeader = TRUE,
        DT::DTOutput(ns("table")),
        shiny::downloadButton(ns("table_download"), "Download prevalence table")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Shared and unique core taxa", width = 12, status = "primary", solidHeader = TRUE,
        shiny::helpText("For 2–4 groups the app draws a Venn diagram. With more groups it draws an UpSet plot."),
        shiny::plotOutput(ns("shared_plot"), height = "650px"),
        shiny::downloadButton(ns("shared_table_download"), "Download group/taxon membership")
      )
    )
  )
}

zamp_taxa_overview_data <- function(ps, rank, group, detection_pct = 0.1) {
  x <- zamp_aggregate_rank(ps, rank)
  mat <- zamp_otu_matrix(x, samples_in_rows = TRUE)
  rel <- zamp_relative_matrix(mat)
  labels <- zamp_rank_labels(x, rank)
  colnames(rel) <- labels
  meta <- zamp_metadata(x)
  rownames(meta) <- meta$Sample
  meta <- meta[rownames(rel), , drop = FALSE]

  groups <- if (identical(group, "None")) factor(rep("All samples", nrow(rel))) else factor(meta[[group]])
  out <- lapply(levels(groups), function(g) {
    idx <- which(groups == g)
    if (!length(idx)) return(NULL)
    r <- rel[idx, , drop = FALSE]
    data.frame(
      Taxon = colnames(r),
      Group = g,
      Prevalence_pct = colMeans(r > detection_pct / 100) * 100,
      Mean_relative_abundance_pct = colMeans(r) * 100,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

zamp_core_sets <- function(ps, rank, group, detection_pct, prevalence_pct) {
  if (identical(group, "None")) return(list())
  x <- zamp_aggregate_rank(ps, rank)
  mat <- zamp_otu_matrix(x, samples_in_rows = TRUE)
  rel <- zamp_relative_matrix(mat)
  labels <- zamp_rank_labels(x, rank)
  colnames(rel) <- labels
  meta <- zamp_metadata(x)
  rownames(meta) <- meta$Sample
  meta <- meta[rownames(rel), , drop = FALSE]
  grp <- factor(meta[[group]])
  sets <- lapply(levels(grp), function(g) {
    r <- rel[grp == g & !is.na(grp), , drop = FALSE]
    if (!nrow(r)) return(character())
    prevalence <- colMeans(r > detection_pct / 100) * 100
    names(prevalence)[prevalence >= prevalence_pct]
  })
  names(sets) <- levels(grp)
  sets[vapply(sets, length, integer(1)) > 0]
}

mod_taxa_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps()))
      vars <- zamp_categorical_variables(ps())
      default_rank <- if ("Genus" %in% ranks) "Genus" else ranks[1]
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = default_rank)
      shiny::updateSelectInput(session, "group", choices = c("None", vars), selected = "None")
    })

    taxa_data <- shiny::eventReactive(input$generate, {
      shiny::req(ps(), input$rank)
      zamp_taxa_overview_data(ps(), input$rank, input$group, input$detection)
    }, ignoreInit = TRUE)

    taxa_plot <- shiny::reactive({
      df <- taxa_data()
      df$Group <- factor(df$Group)
      pal <- stats::setNames(zamp_palette(nlevels(df$Group)), levels(df$Group))
      ggplot2::ggplot(df, ggplot2::aes(
        x = Prevalence_pct,
        y = Mean_relative_abundance_pct,
        colour = Group,
        text = paste0("Taxon: ", Taxon, "<br>Group: ", Group, "<br>Prevalence: ", round(Prevalence_pct, 1), "%<br>Mean abundance: ", signif(Mean_relative_abundance_pct, 3), "%")
      )) +
        ggplot2::geom_point(size = 2.8, alpha = 0.78) +
        ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
        ggplot2::scale_y_continuous(trans = "log1p", labels = scales::label_number(suffix = "%", accuracy = 0.01)) +
        ggplot2::scale_x_continuous(limits = c(0, 100), labels = scales::label_number(suffix = "%")) +
        ggplot2::labs(title = "Taxon prevalence versus abundance", subtitle = paste(input$rank, "level; detection threshold", input$detection, "%"), x = "Prevalence", y = "Mean relative abundance", colour = if (identical(input$group, "None")) NULL else input$group) +
        zamp_theme()
    })

    core_sets <- shiny::reactive({
      shiny::req(taxa_data())
      zamp_core_sets(ps(), input$rank, input$group, input$detection, input$prevalence)
    })

    output$plot <- plotly::renderPlotly({ shiny::req(taxa_plot()); zamp_plotly(taxa_plot(), tooltip = "text") })
    output$table <- DT::renderDT({
      x <- taxa_data()
      x$Prevalence_pct <- round(x$Prevalence_pct, 2)
      x$Mean_relative_abundance_pct <- signif(x$Mean_relative_abundance_pct, 4)
      DT::datatable(x, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE, filter = "top")
    })

    output$shared_plot <- shiny::renderPlot({
      shiny::req(taxa_data())
      shiny::validate(shiny::need(!identical(input$group, "None"), "Choose a grouping variable to calculate shared core taxa."))
      sets <- core_sets()
      shiny::validate(shiny::need(length(sets) >= 2, "At least two groups need core taxa at the selected thresholds."))
      if (length(sets) <= 4) {
        print(
          ggvenn::ggvenn(sets, fill_color = zamp_palette(length(sets)), stroke_size = 0.5, set_name_size = 4) +
            ggplot2::labs(title = paste("Core", input$rank, "shared across", input$group)) + zamp_theme()
        )
      } else {
        upset_df <- UpSetR::fromList(sets)
        UpSetR::upset(
          upset_df,
          nsets = min(length(sets), 12),
          nintersects = 30,
          order.by = "freq",
          keep.order = TRUE,
          sets.bar.color = "#0072B2",
          main.bar.color = "#333333",
          text.scale = 1.15
        )
      }
    })

    zamp_plot_download_server("download", taxa_plot, function() paste0("taxa_prevalence_abundance_", input$rank, "_", Sys.Date()))

    output$table_download <- shiny::downloadHandler(
      filename = function() paste0("taxa_prevalence_", input$rank, "_", Sys.Date(), ".tsv"),
      content = function(file) utils::write.table(taxa_data(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    )
    output$shared_table_download <- shiny::downloadHandler(
      filename = function() paste0("core_taxa_membership_", input$rank, "_", Sys.Date(), ".tsv"),
      content = function(file) {
        sets <- core_sets()
        x <- dplyr::bind_rows(lapply(names(sets), function(g) data.frame(Group = g, Taxon = sets[[g]], stringsAsFactors = FALSE)))
        utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    shiny::observeEvent(input$generate, {
      if (isTRUE(auto_save())) {
        shiny::req(taxa_plot())
        pfile <- zamp_save_plot_to_dir(taxa_plot(), output_dir(), paste0("taxa_prevalence_abundance_", input$rank), "pdf", 7, 5.5)
        tfile <- zamp_save_table_to_dir(taxa_data(), output_dir(), paste0("taxa_prevalence_", input$rank))
        shiny::showNotification(paste("Saved:", basename(pfile), "and", basename(tfile)), type = "message")
      }
    })
  })
}
