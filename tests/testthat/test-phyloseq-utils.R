test_that("phyloseq validation and orientation are stable", {
  ps <- make_test_phyloseq(); expect_silent(zAMPExplorer:::zamp_validate_phyloseq(ps))
  samples_rows <- zAMPExplorer:::zamp_otu_matrix(ps, samples_in_rows = TRUE)
  expect_equal(dim(samples_rows), c(6, 4)); expect_equal(rownames(samples_rows), phyloseq::sample_names(ps))
  taxa_rows <- zAMPExplorer:::zamp_otu_matrix(ps, samples_in_rows = FALSE); expect_equal(dim(taxa_rows), c(4, 6))
})

test_that("transformations return finite sample-by-feature matrices", {
  ps <- make_test_phyloseq(); mat <- zAMPExplorer:::zamp_otu_matrix(ps, samples_in_rows = TRUE)
  rel <- zAMPExplorer:::zamp_transform_matrix(mat, "compositional"); expect_equal(rowSums(rel), rep(1, nrow(rel)))
  clr <- zAMPExplorer:::zamp_transform_matrix(mat, "clr"); expect_true(all(is.finite(clr))); expect_equal(rowMeans(clr), rep(0, nrow(clr)), tolerance = 1e-10)
})

test_that("taxonomic aggregation keeps sample totals", {
  ps <- make_test_phyloseq(); genus <- zAMPExplorer:::zamp_aggregate_rank(ps, "Genus"); expect_equal(phyloseq::sample_sums(genus), phyloseq::sample_sums(ps))
})
