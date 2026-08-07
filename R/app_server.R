#' zAMPExplorer server
#'
#' @param input,output,session Standard Shiny server objects.
#' @noRd
app_server <- function(input, output, session) {
  uploaded_ps <- shiny::reactive({
    shiny::req(input$physeq_file)
    shiny::withProgress(message = "Loading phyloseq object", value = 0.2, {
      ps <- zamp_read_phyloseq(input$physeq_file$datapath)
      shiny::incProgress(0.8)
      ps
    })
  })

  output_dir <- shiny::reactive({
    path <- trimws(input$output_dir %||% "")
    shiny::validate(shiny::need(nzchar(path), "Choose a results folder."))
    path <- path.expand(path)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    shiny::validate(shiny::need(dir.exists(path), paste("Cannot create results folder:", path)))
    normalizePath(path, winslash = "/", mustWork = FALSE)
  })

  auto_save <- shiny::reactive(isTRUE(input$auto_save))

  output$output_dir_status <- shiny::renderText({
    paste("Results folder:", output_dir())
  })

  output$upload_status <- shiny::renderUI({
    shiny::req(uploaded_ps())
    ps <- uploaded_ps()
    components <- c(
      "OTU/ASV table",
      if (tryCatch(!is.null(phyloseq::tax_table(ps)), error = function(e) FALSE)) "taxonomy",
      if (tryCatch(!is.null(phyloseq::sample_data(ps)), error = function(e) FALSE)) "metadata",
      if (tryCatch(!is.null(phyloseq::phy_tree(ps)), error = function(e) FALSE)) "phylogenetic tree",
      if (tryCatch(!is.null(phyloseq::refseq(ps)), error = function(e) FALSE)) "reference sequences"
    )
    shiny::div(
      class = "alert alert-success zamp-upload-status",
      shiny::tags$strong("Phyloseq object ready."),
      shiny::tags$br(),
      sprintf(
        "%s samples · %s taxa · %s reads",
        format(phyloseq::nsamples(ps), big.mark = ","),
        format(phyloseq::ntaxa(ps), big.mark = ","),
        format(sum(phyloseq::sample_sums(ps)), big.mark = ",")
      ),
      shiny::tags$br(),
      paste("Components:", paste(components, collapse = ", "))
    )
  })

  shiny::observeEvent(uploaded_ps(), {
    if (isTRUE(auto_save())) {
      dir <- output_dir()
      manifest <- data.frame(
        item = c("source_filename", "samples", "taxa", "total_reads", "loaded_at"),
        value = c(
          input$physeq_file$name,
          phyloseq::nsamples(uploaded_ps()),
          phyloseq::ntaxa(uploaded_ps()),
          sum(phyloseq::sample_sums(uploaded_ps())),
          base::format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
        ),
        stringsAsFactors = FALSE
      )
      utils::write.table(
        manifest,
        file.path(dir, "zAMPExplorer_session_manifest.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  })

  analysis_ps <- mod_qc_server("qc", uploaded_ps, output_dir, auto_save)
  mod_overview_server("overview", analysis_ps)
  mod_taxa_server("taxa", analysis_ps, output_dir, auto_save)
  mod_composition_server("composition", analysis_ps, output_dir, auto_save)
  mod_heatmap_server("heatmap", analysis_ps, output_dir, auto_save)
  mod_alpha_server("alpha", analysis_ps, output_dir, auto_save)
  mod_beta_server("beta", analysis_ps, output_dir, auto_save)
  mod_maaslin2_server("maaslin2", analysis_ps, output_dir, auto_save)
  mod_dmm_server("dmm", analysis_ps, output_dir, auto_save)
  mod_rda_server("rda", analysis_ps, output_dir, auto_save)
}
