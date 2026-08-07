# Install dependencies needed to run zAMPExplorer from source.

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_cran <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, dependencies = TRUE)
  invisible(missing)
}

install_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
  invisible(missing)
}

cran_packages <- c("DT", "dplyr", "ggplot2", "ggrepel", "ggvenn", "plotly", "scales", "shiny", "shinydashboard", "tidyr", "UpSetR", "writexl")
bioc_packages <- c("ComplexHeatmap", "DirichletMultinomial", "Maaslin2", "phyloseq")
install_cran(cran_packages); install_bioc(bioc_packages); install_cran("vegan")
message("zAMPExplorer dependencies are installed.")
