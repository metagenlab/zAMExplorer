# Phyloseq validation and transformation helpers -----------------------------

zamp_validate_phyloseq <- function(x) {
  if (!inherits(x, "phyloseq")) stop("The uploaded RDS file is not a phyloseq object.")
  otu <- tryCatch(phyloseq::otu_table(x), error = function(e) NULL)
  tax <- tryCatch(phyloseq::tax_table(x), error = function(e) NULL)
  meta <- tryCatch(phyloseq::sample_data(x), error = function(e) NULL)
  if (is.null(otu)) stop("The phyloseq object has no OTU/ASV table.")
  if (is.null(tax)) stop("The phyloseq object has no taxonomy table.")
  if (is.null(meta)) stop("The phyloseq object has no sample metadata.")
  invisible(TRUE)
}

zamp_prepare_phyloseq <- function(x) {
  zamp_validate_phyloseq(x)
  x <- phyloseq::prune_samples(phyloseq::sample_sums(x) > 0, x)
  x <- phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x)
  x
}

zamp_read_phyloseq <- function(path) {
  x <- readRDS(path)
  zamp_prepare_phyloseq(x)
}

zamp_metadata <- function(ps) {
  out <- as(phyloseq::sample_data(ps), "data.frame")
  out$Sample <- rownames(out)
  out
}

zamp_otu_matrix <- function(ps, samples_in_rows = TRUE) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps) == samples_in_rows) mat <- t(mat)
  storage.mode(mat) <- "double"
  mat
}

zamp_aggregate_rank <- function(ps, rank) {
  if (identical(rank, "ASV")) return(ps)
  if (!rank %in% phyloseq::rank_names(ps)) stop("Unknown taxonomic rank: ", rank)
  phyloseq::tax_glom(ps, taxrank = rank, NArm = FALSE)
}

zamp_rank_labels <- function(ps, rank) {
  if (identical(rank, "ASV")) return(phyloseq::taxa_names(ps))
  tt <- as(phyloseq::tax_table(ps), "matrix")
  labels <- as.character(tt[, rank])
  missing <- is.na(labels) | !nzchar(labels)
  labels[missing] <- paste0("Unclassified_", rank, "_", seq_len(sum(missing)))
  labels
}

zamp_taxon_labels <- function(ps, rank) {
  data.frame(
    TaxonID = phyloseq::taxa_names(ps),
    Taxon = zamp_rank_labels(ps, rank),
    stringsAsFactors = FALSE
  )
}

zamp_relative_matrix <- function(mat) {
  rs <- rowSums(mat)
  rs[rs == 0] <- 1
  mat / rs
}

zamp_transform_matrix <- function(mat, method = "compositional", pseudocount = 0.5) {
  method <- tolower(method)
  if (method == "identity") return(mat)
  if (method == "compositional") return(zamp_relative_matrix(mat))
  if (method == "binary") return((mat > 0) * 1)
  if (method == "hellinger") return(sqrt(zamp_relative_matrix(mat)))
  if (method == "log10p") return(log10(mat + 1))
  if (method == "clr") {
    x <- mat + pseudocount
    lx <- log(x)
    return(lx - rowMeans(lx))
  }
  if (method == "z") {
    z <- scale(mat)
    z[!is.finite(z)] <- 0
    return(z)
  }
  stop("Unsupported transformation: ", method)
}

zamp_distance <- function(ps, rank, metric = "bray", transform = "compositional") {
  ps_rank <- zamp_aggregate_rank(ps, rank)
  mat <- zamp_otu_matrix(ps_rank, samples_in_rows = TRUE)
  metric <- tolower(metric)

  if (metric == "aitchison") {
    x <- zamp_transform_matrix(mat, "clr")
    return(stats::dist(x, method = "euclidean"))
  }
  if (metric == "jaccard") {
    return(vegan::vegdist(mat, method = "jaccard", binary = TRUE))
  }
  if (metric == "bray") {
    x <- zamp_transform_matrix(mat, transform)
    return(vegan::vegdist(x, method = "bray"))
  }
  stop("Unsupported distance metric: ", metric)
}

zamp_categorical_variables <- function(ps) {
  meta <- as(phyloseq::sample_data(ps), "data.frame")
  names(meta)[vapply(meta, function(x) is.factor(x) || is.character(x) || is.logical(x), logical(1))]
}

zamp_all_metadata_variables <- function(ps) {
  colnames(as(phyloseq::sample_data(ps), "data.frame"))
}

zamp_complete_metadata <- function(ps, variables) {
  meta <- as(phyloseq::sample_data(ps), "data.frame")
  variables <- intersect(variables, colnames(meta))
  if (!length(variables)) return(meta[FALSE, , drop = FALSE])
  meta[stats::complete.cases(meta[, variables, drop = FALSE]), , drop = FALSE]
}
