mod_overview_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboard::valueBoxOutput(ns("samples"), width = 3),
      shinydashboard::valueBoxOutput(ns("taxa"), width = 3),
      shinydashboard::valueBoxOutput(ns("reads"), width = 3),
      shinydashboard::valueBoxOutput(ns("variables"), width = 3)
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Sample metadata", width = 12, status = "primary", solidHeader = TRUE,
        DT::DTOutput(ns("metadata")),
        shiny::downloadButton(ns("metadata_download"), "Download metadata")
      )
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Abundance + taxonomy", width = 12, status = "primary", solidHeader = TRUE,
        DT::DTOutput(ns("combined")),
        shiny::selectInput(ns("table_format"), "Download format", c("TSV" = "tsv", "CSV" = "csv", "Excel" = "xlsx")),
        shiny::downloadButton(ns("combined_download"), "Download table")
      )
    )
  )
}

mod_overview_server <- function(id, ps) {
  shiny::moduleServer(id, function(input, output, session) {
    summary_vals <- shiny::reactive({
      shiny::req(ps())
      c(
        samples = phyloseq::nsamples(ps()),
        taxa = phyloseq::ntaxa(ps()),
        reads = sum(phyloseq::sample_sums(ps())),
        variables = length(phyloseq::sample_variables(ps()))
      )
    })

    output$samples <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(format(summary_vals()[["samples"]], big.mark = ","), "Samples", icon = shiny::icon("vials"), color = "aqua")
    })
    output$taxa <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(format(summary_vals()[["taxa"]], big.mark = ","), "ASVs / taxa", icon = shiny::icon("dna"), color = "green")
    })
    output$reads <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(format(summary_vals()[["reads"]], big.mark = ","), "Total reads", icon = shiny::icon("layer-group"), color = "yellow")
    })
    output$variables <- shinydashboard::renderValueBox({
      shinydashboard::valueBox(summary_vals()[["variables"]], "Metadata variables", icon = shiny::icon("table"), color = "purple")
    })

    metadata_df <- shiny::reactive({
      shiny::req(ps())
      zamp_metadata(ps())
    })

    combined_df <- shiny::reactive({
      shiny::req(ps())
      otu <- zamp_otu_matrix(ps(), samples_in_rows = FALSE)
      otu <- data.frame(ASV = rownames(otu), otu, check.names = FALSE)
      tax <- as.data.frame(phyloseq::tax_table(ps()), stringsAsFactors = FALSE)
      tax$ASV <- rownames(tax)
      dplyr::left_join(tax, otu, by = "ASV")
    })

    output$metadata <- DT::renderDT({
      DT::datatable(metadata_df(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE, filter = "top")
    })
    output$combined <- DT::renderDT({
      DT::datatable(combined_df(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE, filter = "top")
    })

    output$metadata_download <- shiny::downloadHandler(
      filename = function() paste0("zAMPExplorer_metadata_", Sys.Date(), ".tsv"),
      content = function(file) utils::write.table(metadata_df(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    )
    output$combined_download <- shiny::downloadHandler(
      filename = function() paste0("zAMPExplorer_abundance_taxonomy_", Sys.Date(), ".", input$table_format),
      content = function(file) {
        x <- combined_df()
        if (input$table_format == "xlsx") writexl::write_xlsx(x, file)
        else if (input$table_format == "csv") utils::write.csv(x, file, row.names = FALSE)
        else utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  })
}
