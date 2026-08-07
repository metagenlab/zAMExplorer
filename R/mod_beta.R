mod_beta_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::box(title = "Ordination settings", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("method"), "Ordination", c("PCoA", "PCA", "NMDS")),
        shiny::selectInput(ns("distance"), "Distance", c("Bray-Curtis" = "bray", "Jaccard" = "jaccard", "Aitchison" = "aitchison")),
        shiny::selectInput(ns("transform"), "Transformation", c("Compositional" = "compositional", "Hellinger" = "hellinger", "CLR" = "clr", "Identity" = "identity", "Binary" = "binary", "log10(x+1)" = "log10p")),
        shiny::selectInput(ns("rank"), "Taxonomic rank", choices = NULL),
        shiny::selectInput(ns("colour"), "Colour by", choices = NULL),
        shiny::selectInput(ns("shape"), "Shape by", choices = NULL),
        shiny::checkboxInput(ns("ellipse"), "Show 95% group ellipses", TRUE),
        shiny::actionButton(ns("generate"), "Generate ordination", class = "btn-primary")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Beta-diversity ordination", width = 12, status = "primary", solidHeader = TRUE,
        plotly::plotlyOutput(ns("plot"), height = "680px"),
        shiny::verbatimTextOutput(ns("note")),
        zamp_plot_download_ui(ns("download"), 7, 6)
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "PERMANOVA and dispersion", width = 12, status = "primary", solidHeader = TRUE,
        shiny::selectInput(ns("stats_group"), "Grouping variable", choices = NULL),
        shiny::selectInput(ns("strata"), "Permutation strata", choices = NULL),
        shiny::numericInput(ns("permutations"), "Permutations", 999, min = 99, max = 99999, step = 100),
        shiny::actionButton(ns("run_stats"), "Run PERMANOVA + dispersion", class = "btn-primary"),
        shiny::h4("PERMANOVA"), DT::DTOutput(ns("permanova")),
        shiny::h4("Dispersion"), DT::DTOutput(ns("dispersion")),
        shiny::downloadButton(ns("stats_download"), "Download beta-diversity statistics")
      )
    )
  )
}

zamp_ordination <- function(ps, rank, method, distance, transform, colour_var, shape_var = "None", ellipse = TRUE) {
  ps_rank <- zamp_aggregate_rank(ps, rank)
  mat <- zamp_otu_matrix(ps_rank, samples_in_rows = TRUE)
  meta <- zamp_metadata(ps_rank); rownames(meta) <- meta$Sample; meta <- meta[rownames(mat), , drop = FALSE]
  if (distance == "aitchison") transform <- "clr"
  if (distance == "jaccard") transform <- "binary"
  x <- zamp_transform_matrix(mat, transform)
  note <- NULL
  if (method == "PCA") {
    fit <- stats::prcomp(x, center = TRUE, scale. = FALSE)
    coords <- fit$x[, 1:2, drop = FALSE]
    ve <- fit$sdev^2 / sum(fit$sdev^2) * 100
    xlab <- sprintf("PC1 (%.1f%%)", ve[1]); ylab <- sprintf("PC2 (%.1f%%)", ve[2])
  } else {
    dist_obj <- if (distance == "aitchison") stats::dist(x) else if (distance == "jaccard") vegan::vegdist(mat, "jaccard", binary = TRUE) else vegan::vegdist(x, "bray")
    if (method == "PCoA") {
      fit <- stats::cmdscale(dist_obj, k = 2, eig = TRUE, add = TRUE)
      coords <- fit$points[, 1:2, drop = FALSE]
      pos <- fit$eig[fit$eig > 0]; ve <- if (length(pos)) fit$eig[1:2] / sum(pos) * 100 else c(NA_real_, NA_real_)
      xlab <- if (is.finite(ve[1])) sprintf("PCoA1 (%.1f%%)", ve[1]) else "PCoA1"
      ylab <- if (is.finite(ve[2])) sprintf("PCoA2 (%.1f%%)", ve[2]) else "PCoA2"
    } else {
      fit <- vegan::metaMDS(dist_obj, k = 2, trymax = 80, autotransform = FALSE, trace = FALSE)
      coords <- vegan::scores(fit, display = "sites", choices = 1:2)
      xlab <- "NMDS1"; ylab <- "NMDS2"; note <- sprintf("NMDS stress = %.3f", fit$stress)
    }
  }
  df <- data.frame(Sample = rownames(coords), Axis1 = coords[, 1], Axis2 = coords[, 2], stringsAsFactors = FALSE)
  df <- dplyr::left_join(df, meta, by = "Sample"); df$Colour <- factor(df[[colour_var]])
  if (!identical(shape_var, "None")) df$Shape <- factor(df[[shape_var]])
  pal <- stats::setNames(zamp_palette(nlevels(df$Colour)), levels(df$Colour))
  aes_base <- if (identical(shape_var, "None")) ggplot2::aes(x = Axis1, y = Axis2, colour = Colour, text = Sample) else ggplot2::aes(x = Axis1, y = Axis2, colour = Colour, shape = Shape, text = Sample)
  p <- ggplot2::ggplot(df, aes_base) + ggplot2::geom_point(size = 3.2, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
    ggplot2::labs(title = paste(method, "ordination"), subtitle = paste(rank, "level ·", distance, "distance ·", transform), x = xlab, y = ylab, colour = colour_var, shape = if (identical(shape_var, "None")) NULL else shape_var) +
    zamp_theme() + ggplot2::theme(aspect.ratio = 1)
  if (isTRUE(ellipse)) {
    counts <- table(df$Colour); valid <- names(counts[counts >= 3])
    if (length(valid)) p <- p + ggplot2::stat_ellipse(data = df[df$Colour %in% valid, , drop = FALSE], ggplot2::aes(group = Colour, colour = Colour), level = 0.95, linewidth = 0.65, alpha = 0.8, show.legend = FALSE)
  }
  list(plot = p, data = df, note = note)
}

mod_beta_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps())); vars <- zamp_categorical_variables(ps())
      default_rank <- if ("Genus" %in% ranks) "Genus" else ranks[1]
      shiny::updateSelectInput(session, "rank", choices = ranks, selected = default_rank)
      shiny::updateSelectInput(session, "colour", choices = vars, selected = vars[1])
      shiny::updateSelectInput(session, "shape", choices = c("None", vars), selected = "None")
      shiny::updateSelectInput(session, "stats_group", choices = vars, selected = vars[1])
      shiny::updateSelectInput(session, "strata", choices = c("None", vars), selected = "None")
    })
    ord <- shiny::eventReactive(input$generate, { shiny::req(ps(), input$rank, input$colour); zamp_ordination(ps(), input$rank, input$method, input$distance, input$transform, input$colour, input$shape, input$ellipse) }, ignoreInit = TRUE)
    ord_plot <- shiny::reactive({ shiny::req(ord()); ord()$plot })
    output$plot <- plotly::renderPlotly({ shiny::req(ord_plot()); zamp_plotly(ord_plot(), tooltip = c("text", "x", "y", "colour")) })
    output$note <- shiny::renderText({ shiny::req(ord()); ord()$note })
    zamp_plot_download_server("download", ord_plot, function() paste0("beta_", input$method, "_", input$distance, "_", Sys.Date()))

    beta_stats <- shiny::eventReactive(input$run_stats, {
      shiny::req(ps(), input$rank, input$stats_group)
      dist_obj <- zamp_distance(ps(), input$rank, input$distance, input$transform)
      meta <- zamp_metadata(ps()); rownames(meta) <- meta$Sample; ids <- attr(dist_obj, "Labels"); meta <- meta[ids, , drop = FALSE]
      group <- factor(meta[[input$stats_group]]); keep <- !is.na(group); meta <- meta[keep, , drop = FALSE]; group <- droplevels(group[keep])
      mat_dist <- as.matrix(dist_obj)[keep, keep, drop = FALSE]; dist_use <- stats::as.dist(mat_dist)
      shiny::validate(shiny::need(nlevels(group) >= 2, "At least two groups are required."))
      strata <- NULL; if (!identical(input$strata, "None")) strata <- meta[[input$strata]]
      dat <- data.frame(group = group)
      perm <- vegan::adonis2(dist_use ~ group, data = dat, permutations = input$permutations, strata = strata)
      perm_df <- data.frame(Term = rownames(perm), as.data.frame(perm), row.names = NULL, check.names = FALSE)
      bd <- vegan::betadisper(dist_use, group); bd_anova <- stats::anova(bd); bd_perm <- vegan::permutest(bd, permutations = input$permutations)
      dispersion_df <- data.frame(Test = c("ANOVA", "Permutation"), F = c(unname(bd_anova$`F value`[1]), unname(bd_perm$tab[1, "F"])), p_value = c(unname(bd_anova$`Pr(>F)`[1]), unname(bd_perm$tab[1, "Pr(>F)"])))
      list(permanova = perm_df, dispersion = dispersion_df)
    }, ignoreInit = TRUE)
    output$permanova <- DT::renderDT({ shiny::req(beta_stats()); DT::datatable(beta_stats()$permanova, options = list(dom = "t"), rownames = FALSE) })
    output$dispersion <- DT::renderDT({ shiny::req(beta_stats()); DT::datatable(beta_stats()$dispersion, options = list(dom = "t"), rownames = FALSE) })
    output$stats_download <- shiny::downloadHandler(filename = function() paste0("beta_diversity_stats_", Sys.Date(), ".tsv"), content = function(file) {
      shiny::req(beta_stats()); x <- dplyr::bind_rows(cbind(Analysis = "PERMANOVA", beta_stats()$permanova), cbind(Analysis = "Dispersion", beta_stats()$dispersion)); utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
    })
    shiny::observeEvent(input$generate, { if (isTRUE(auto_save())) { shiny::req(ord_plot()); file <- zamp_save_plot_to_dir(ord_plot(), output_dir(), paste0("beta_", input$method, "_", input$distance), "pdf", 7, 6); shiny::showNotification(paste("Saved:", file), type = "message") } })
    shiny::observeEvent(input$run_stats, { if (isTRUE(auto_save())) { shiny::req(beta_stats()); zamp_save_table_to_dir(beta_stats()$permanova, output_dir(), "beta_PERMANOVA"); zamp_save_table_to_dir(beta_stats()$dispersion, output_dir(), "beta_dispersion") } })
  })
}
