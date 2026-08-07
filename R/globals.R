`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

utils::globalVariables(c(
  ".data",
  "Association",
  "Axis1",
  "Axis2",
  "Colour",
  "Component",
  "Criterion",
  "Depth",
  "Diversity",
  "FeatureAssociation",
  "Mean_relative_abundance_pct",
  "Prevalence_pct",
  "Richness",
  "Shape",
  "Significant",
  "Taxon",
  "TaxonID",
  "Taxon_display",
  "Value",
  "Weight",
  "coef",
  "feature",
  "k"
))

#' @importFrom methods as
#' @importFrom stats coef
NULL