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

  shiny::validate(
    shiny::need(nrow(mat) >= 3, "At least three samples are required for a two-dimensional ordination."),
    shiny::need(ncol(mat) >= 2, "At least two taxa/features are required for ordination.")
  )

  meta <- zamp_metadata(ps_rank)
  rownames(meta) <- meta$Sample
  meta <- meta[rownames(mat), , drop = FALSE]

  shiny::validate(
    shiny::need(colour_var %in% colnames(meta), "The selected colour variable is not available in the metadata."),
    shiny::need(identical(shape_var, "None") || shape_var %in% colnames(meta), "The selected shape variable is not available in the metadata.")
  )

  if (distance == "aitchison") transform <- "clr"
  if (distance == "jaccard") transform <- "binary"
  x <- zamp_transform_matrix(mat, transform)

  shiny::validate(
    shiny::need(all(is.finite(x)), "The selected transformation produced non-finite values.")
  )

  note <- NULL

  if (method == "PCA") {
    informative <- apply(x, 2, function(v) {
      s <- stats::sd(v, na.rm = TRUE)
      is.finite(s) && s > 0
    })
    x_pca <- x[, informative, drop = FALSE]

    shiny::validate(
      shiny::need(ncol(x_pca) >= 2, "PCA requires at least two taxa/features with non-zero variance.")
    )

    fit <- stats::prcomp(x_pca, center = TRUE, scale. = FALSE)
    shiny::validate(
      shiny::need(ncol(fit$x) >= 2, "PCA produced fewer than two ordination axes.")
    )

    coords <- fit$x[, 1:2, drop = FALSE]
    ve <- fit$sdev^2 / sum(fit$sdev^2) * 100
    xlab <- sprintf("PC1 (%.1f%%)", ve[1])
    ylab <- sprintf("PC2 (%.1f%%)", ve[2])
  } else {
    dist_obj <- if (distance == "aitchison") {
      stats::dist(x)
    } else if (distance == "jaccard") {
      vegan::vegdist(mat, "jaccard", binary = TRUE)
    } else {
      vegan::vegdist(x, "bray")
    }

    dist_values <- as.numeric(dist_obj)
    shiny::validate(
      shiny::need(length(dist_values) > 0 && all(is.finite(dist_values)), "Distance calculation returned non-finite values."),
      shiny::need(any(dist_values > 0), "All samples have zero pairwise distance; ordination cannot be calculated.")
    )

    if (method == "PCoA") {
      fit <- stats::cmdscale(dist_obj, k = 2, eig = TRUE, add = TRUE)
      shiny::validate(
        shiny::need(ncol(as.matrix(fit$points)) >= 2, "PCoA produced fewer than two ordination axes.")
      )

      coords <- fit$points[, 1:2, drop = FALSE]
      pos <- fit$eig[fit$eig > 0]
      ve <- if (length(pos)) fit$eig[1:2] / sum(pos) * 100 else c(NA_real_, NA_real_)
      xlab <- if (is.finite(ve[1])) sprintf("PCoA1 (%.1f%%)", ve[1]) else "PCoA1"
      ylab <- if (is.finite(ve[2])) sprintf("PCoA2 (%.1f%%)", ve[2]) else "PCoA2"
    } else {
      fit <- vegan::metaMDS(
        dist_obj,
        k = 2,
        trymax = 80,
        autotransform = FALSE,
        trace = FALSE
      )
      coords <- vegan::scores(fit, display = "sites", choices = 1:2)
      shiny::validate(
        shiny::need(ncol(as.matrix(coords)) >= 2, "NMDS produced fewer than two ordination axes.")
      )

      xlab <- "NMDS1"
      ylab <- "NMDS2"
      note <- sprintf("NMDS stress = %.3f", fit$stress)
    }
  }

  df <- data.frame(
    Sample = rownames(coords),
    Axis1 = coords[, 1],
    Axis2 = coords[, 2],
    stringsAsFactors = FALSE
  )
  df <- dplyr::left_join(df, meta, by = "Sample")
  df$Colour <- factor(df[[colour_var]])
  df <- df[!is.na(df$Colour), , drop = FALSE]
  df$Colour <- droplevels(df$Colour)

  shiny::validate(
    shiny::need(nrow(df) >= 2, "Too few samples remain after removing missing values in the colour variable."),
    shiny::need(nlevels(df$Colour) >= 1, "No non-missing colour groups remain.")
  )

  if (!identical(shape_var, "None")) {
    df$Shape <- factor(df[[shape_var]])
    df <- df[!is.na(df$Shape), , drop = FALSE]
    df$Shape <- droplevels(df$Shape)
    shiny::validate(
      shiny::need(nrow(df) >= 2, "Too few samples remain after removing missing values in the shape variable."),
      shiny::need(
        nlevels(df$Shape) <= 12,
        paste0(
          "The selected shape variable has ", nlevels(df$Shape),
          " levels. Choose 'None' or a variable with 12 or fewer levels for a readable ordination."
        )
      )
    )
  }

  pal <- stats::setNames(zamp_palette(nlevels(df$Colour)), levels(df$Colour))
  aes_base <- if (identical(shape_var, "None")) {
    ggplot2::aes(x = Axis1, y = Axis2, colour = Colour, text = Sample)
  } else {
    ggplot2::aes(x = Axis1, y = Axis2, colour = Colour, shape = Shape, text = Sample)
  }

  p <- ggplot2::ggplot(df, aes_base) +
    ggplot2::geom_point(size = 3.2, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE)

  if (!identical(shape_var, "None")) {
    shape_values <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6, 9)
    p <- p + ggplot2::scale_shape_manual(values = shape_values[seq_len(nlevels(df$Shape))], drop = FALSE)
  }

  p <- p +
    ggplot2::labs(
      title = paste(method, "ordination"),
      subtitle = paste(rank, "level ·", distance, "distance ·", transform),
      x = xlab,
      y = ylab,
      colour = colour_var,
      shape = if (identical(shape_var, "None")) NULL else shape_var
    ) +
    zamp_theme() +
    ggplot2::theme(aspect.ratio = 1)

  if (isTRUE(ellipse)) {
    counts <- table(df$Colour)
    valid <- names(counts[counts >= 3])
    if (length(valid)) {
      p <- p + ggplot2::stat_ellipse(
        data = df[df$Colour %in% valid, , drop = FALSE],
        ggplot2::aes(group = Colour, colour = Colour),
        level = 0.95,
        linewidth = 0.65,
        alpha = 0.8,
        show.legend = FALSE
      )
    }
  }

  list(plot = p, data = df, note = note)
}

mod_beta_server <- function(id, ps, output_dir, auto_save) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(ps(), {
      shiny::req(ps())
      ranks <- c("ASV", phyloseq::rank_names(ps()))
      vars <- zamp_categorical_variables(ps())
      default_rank <- if ("Genus" %in% ranks) "Genus" else ranks[1]

      shiny::updateSelectInput(session, "rank", choices = ranks, selected = default_rank)
      shiny::updateSelectInput(session, "colour", choices = vars, selected = if (length(vars)) vars[1] else character(0))
      shiny::updateSelectInput(session, "shape", choices = c("None", vars), selected = "None")
      shiny::updateSelectInput(session, "stats_group", choices = vars, selected = if (length(vars)) vars[1] else character(0))
      shiny::updateSelectInput(session, "strata", choices = c("None", vars), selected = "None")
    })

    ord <- shiny::eventReactive(input$generate, {
      shiny::req(ps(), input$rank, input$colour)
      tryCatch(
        zamp_ordination(
          ps(), input$rank, input$method, input$distance, input$transform,
          input$colour, input$shape, input$ellipse
        ),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, paste("Beta-diversity ordination failed:", conditionMessage(e))))
        }
      )
    }, ignoreInit = TRUE)

    ord_plot <- shiny::reactive({
      shiny::req(ord())
      ord()$plot
    })

    output$plot <- plotly::renderPlotly({
      shiny::req(ord_plot())
      zamp_plotly(ord_plot(), tooltip = c("text", "x", "y", "colour"))
    })

    output$note <- shiny::renderText({
      shiny::req(ord())
      ord()$note
    })

    zamp_plot_download_server(
      "download",
      ord_plot,
      function() paste0("beta_", input$method, "_", input$distance, "_", Sys.Date())
    )

    beta_stats <- shiny::eventReactive(input$run_stats, {
      shiny::req(ps(), input$rank, input$stats_group)

      tryCatch({
        meta <- zamp_metadata(ps())
        rownames(meta) <- meta$Sample

        vars_for_complete_cases <- input$stats_group
        if (!identical(input$strata, "None")) {
          vars_for_complete_cases <- unique(c(vars_for_complete_cases, input$strata))
        }

        shiny::validate(
          shiny::need(
            all(vars_for_complete_cases %in% colnames(meta)),
            "The selected PERMANOVA grouping/strata variable is not available in the metadata."
          )
        )

        dist_obj <- zamp_distance(ps(), input$rank, input$distance, input$transform)
        ids <- attr(dist_obj, "Labels")
        meta <- meta[ids, , drop = FALSE]

        keep <- stats::complete.cases(meta[, vars_for_complete_cases, drop = FALSE])
        meta <- meta[keep, , drop = FALSE]

        shiny::validate(
          shiny::need(nrow(meta) >= 3, "At least three complete samples are required for PERMANOVA and dispersion.")
        )

        group <- droplevels(factor(meta[[input$stats_group]]))
        shiny::validate(
          shiny::need(nlevels(group) >= 2, "At least two groups are required.")
        )

        mat_dist <- as.matrix(dist_obj)[keep, keep, drop = FALSE]
        dist_use <- stats::as.dist(mat_dist)

        strata <- NULL
        if (!identical(input$strata, "None")) {
          strata <- factor(meta[[input$strata]])
          shiny::validate(
            shiny::need(nlevels(droplevels(strata)) >= 2, "The selected permutation strata has fewer than two non-missing levels.")
          )
        }

        dat <- data.frame(group = group)
        perm <- vegan::adonis2(
          dist_use ~ group,
          data = dat,
          permutations = input$permutations,
          strata = strata
        )
        perm_df <- data.frame(
          Term = rownames(perm),
          as.data.frame(perm),
          row.names = NULL,
          check.names = FALSE
        )

        bd <- vegan::betadisper(dist_use, group)
        bd_anova <- stats::anova(bd)
        bd_perm <- vegan::permutest(bd, permutations = input$permutations)
        dispersion_df <- data.frame(
          Test = c("ANOVA", "Permutation"),
          F = c(
            unname(bd_anova$`F value`[1]),
            unname(bd_perm$tab[1, "F"])
          ),
          p_value = c(
            unname(bd_anova$`Pr(>F)`[1]),
            unname(bd_perm$tab[1, "Pr(>F)"])
          )
        )

        list(permanova = perm_df, dispersion = dispersion_df)
      }, error = function(e) {
        shiny::validate(shiny::need(FALSE, paste("Beta-diversity statistics failed:", conditionMessage(e))))
      })
    }, ignoreInit = TRUE)

    output$permanova <- DT::renderDT({
      shiny::req(beta_stats())
      DT::datatable(beta_stats()$permanova, options = list(dom = "t"), rownames = FALSE)
    })

    output$dispersion <- DT::renderDT({
      shiny::req(beta_stats())
      DT::datatable(beta_stats()$dispersion, options = list(dom = "t"), rownames = FALSE)
    })

    output$stats_download <- shiny::downloadHandler(
      filename = function() paste0("beta_diversity_stats_", Sys.Date(), ".tsv"),
      content = function(file) {
        shiny::req(beta_stats())
        x <- dplyr::bind_rows(
          cbind(Analysis = "PERMANOVA", beta_stats()$permanova),
          cbind(Analysis = "Dispersion", beta_stats()$dispersion)
        )
        utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    shiny::observeEvent(input$generate, {
      if (!isTRUE(auto_save())) return()
      tryCatch({
        shiny::req(ord_plot())
        file <- zamp_save_plot_to_dir(
          ord_plot(), output_dir(),
          paste0("beta_", input$method, "_", input$distance),
          "pdf", 7, 6
        )
        shiny::showNotification(paste("Saved:", basename(file)), type = "message")
      }, error = function(e) {
        shiny::showNotification(
          paste("Beta-diversity plot was generated, but automatic saving failed:", conditionMessage(e)),
          type = "warning",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$run_stats, {
      if (!isTRUE(auto_save())) return()
      tryCatch({
        shiny::req(beta_stats())
        zamp_save_table_to_dir(beta_stats()$permanova, output_dir(), "beta_PERMANOVA")
        zamp_save_table_to_dir(beta_stats()$dispersion, output_dir(), "beta_dispersion")
      }, error = function(e) {
        shiny::showNotification(
          paste("Beta-diversity statistics ran, but automatic saving failed:", conditionMessage(e)),
          type = "warning",
          duration = 8
        )
      })
    }, ignoreInit = TRUE)
  })
}
