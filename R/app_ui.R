#' zAMPExplorer user interface
#'
#' @param request Internal Shiny request object.
#' @noRd
app_ui <- function(request) {
  shinydashboard::dashboardPage(
    skin = "blue",
    shinydashboard::dashboardHeader(
      title = shiny::tagList(shiny::icon("dna"), "zAMPExplorer"),
      titleWidth = 260,
      shiny::tags$li(
        class = "dropdown",
        shiny::tags$a(
          href = "https://github.com/metagenlab/zAMPExplorer",
          target = "_blank",
          rel = "noopener noreferrer",
          shiny::icon("book"),
          " Documentation"
        )
      )
    ),
    shinydashboard::dashboardSidebar(
      width = 260,
      shinydashboard::sidebarMenu(
        id = "main_tabs",
        shinydashboard::menuItem("Upload & session", tabName = "upload", icon = shiny::icon("upload")),
        shinydashboard::menuItem("Dataset overview", tabName = "overview", icon = shiny::icon("table")),
        shinydashboard::menuItem("Read QC", tabName = "qc", icon = shiny::icon("chart-column")),
        shinydashboard::menuItem("Taxa overview", tabName = "taxa", icon = shiny::icon("magnifying-glass-chart")),
        shinydashboard::menuItem("Composition", tabName = "composition", icon = shiny::icon("chart-area")),
        shinydashboard::menuItem("Heatmap", tabName = "heatmap", icon = shiny::icon("th")),
        shinydashboard::menuItem("Alpha diversity", tabName = "alpha", icon = shiny::icon("seedling")),
        shinydashboard::menuItem("Beta diversity", tabName = "beta", icon = shiny::icon("circle-nodes")),
        shinydashboard::menuItem("Differential abundance", tabName = "maaslin2", icon = shiny::icon("sliders")),
        shinydashboard::menuItem("Community typing (DMM)", tabName = "dmm", icon = shiny::icon("layer-group")),
        shinydashboard::menuItem("RDA / dbRDA", tabName = "rda", icon = shiny::icon("arrows-to-dot"))
      )
    ),
    shinydashboard::dashboardBody(
      shiny::tags$head(
        shiny::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
        shiny::includeCSS(app_sys("app", "www", "zamp.css"))
      ),
      shinydashboard::tabItems(
        shinydashboard::tabItem(
          tabName = "upload",
          shiny::fluidRow(
            shinydashboard::box(
              title = "1. Load a phyloseq object",
              width = 7,
              status = "primary",
              solidHeader = TRUE,
              shiny::fileInput(
                "physeq_file",
                "Phyloseq RDS file",
                accept = c(".rds", "application/octet-stream"),
                width = "100%"
              ),
              shiny::helpText(
                "The uploaded object is validated, zero-depth samples and zero-abundance taxa are removed, and the cleaned object is used throughout the app."
              ),
              shiny::uiOutput("upload_status")
            ),
            shinydashboard::box(
              title = "2. Publication output",
              width = 5,
              status = "primary",
              solidHeader = TRUE,
              shiny::textInput(
                "output_dir",
                "Results folder on the machine running zAMPExplorer",
                value = file.path(getwd(), "zAMPExplorer_results"),
                width = "100%"
              ),
              shiny::checkboxInput(
                "auto_save",
                "Automatically save generated publication figures and result tables",
                value = TRUE
              ),
              shiny::helpText(
                "Figures are saved as vector PDF by default; PNG downloads use 600 dpi. The same underlying plot object is used on screen and for export."
              ),
              shiny::verbatimTextOutput("output_dir_status")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Analysis workflow",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              shiny::tags$ol(
                shiny::tags$li("Upload a phyloseq .rds object."),
                shiny::tags$li("Inspect metadata and sequencing depth; optionally filter low-depth samples."),
                shiny::tags$li("Explore taxa composition and diversity."),
                shiny::tags$li("Run beta-diversity statistics, DMM, or constrained ordination as appropriate."),
                shiny::tags$li("Use the PDF/SVG/600-dpi PNG controls to export publication-ready figures.")
              ),
              shiny::tags$p(
                class = "zamp-note",
                "Statistical choices remain the responsibility of the analyst. zAMPExplorer records the selected settings in plot subtitles and exported filenames where practical."
              )
            )
          )
        ),
        shinydashboard::tabItem(tabName = "overview", mod_overview_ui("overview")),
        shinydashboard::tabItem(tabName = "qc", mod_qc_ui("qc")),
        shinydashboard::tabItem(tabName = "taxa", mod_taxa_ui("taxa")),
        shinydashboard::tabItem(tabName = "composition", mod_composition_ui("composition")),
        shinydashboard::tabItem(tabName = "heatmap", mod_heatmap_ui("heatmap")),
        shinydashboard::tabItem(tabName = "alpha", mod_alpha_ui("alpha")),
        shinydashboard::tabItem(tabName = "beta", mod_beta_ui("beta")),
        shinydashboard::tabItem(tabName = "maaslin2", mod_maaslin2_ui("maaslin2")),
        shinydashboard::tabItem(tabName = "dmm", mod_dmm_ui("dmm")),
        shinydashboard::tabItem(tabName = "rda", mod_rda_ui("rda"))
      )
    )
  )
}
