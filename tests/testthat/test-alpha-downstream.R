test_that("Alpha Diversity keeps freshly calculated metrics when metadata already contains them", {
  ps <- make_test_phyloseq()

  meta <- as(phyloseq::sample_data(ps), "data.frame")
  meta$Observed <- 999
  meta$Shannon <- 999
  meta$Chao1 <- 999
  phyloseq::sample_data(ps) <- phyloseq::sample_data(meta)

  res <- zamp_alpha_calculate(
    ps,
    measures = c("Observed", "Shannon", "Chao1")
  )

  expect_true(all(c("Observed", "Shannon", "Chao1") %in% colnames(res$table)))
  expect_false(any(c("Observed.x", "Observed.y", "Shannon.x", "Shannon.y", "Chao1.x", "Chao1.y") %in% colnames(res$table)))
  expect_setequal(
    res$replaced_metadata_metrics,
    c("Observed", "Shannon", "Chao1")
  )

  expected <- suppressWarnings(
    phyloseq::estimate_richness(
      ps,
      measures = c("Observed", "Shannon", "Chao1")
    )
  )

  expect_equal(res$table$Observed, expected[res$table$Sample, "Observed"])
  expect_equal(res$table$Shannon, expected[res$table$Sample, "Shannon"])
  expect_equal(res$table$Chao1, expected[res$table$Sample, "Chao1"])
})

test_that("Beta ordination returns a publication ggplot and records effective Aitchison CLR", {
  ps <- make_test_phyloseq()

  res <- zamp_ordination(
    ps = ps,
    rank = "Genus",
    method = "PCoA",
    distance = "aitchison",
    transform = "hellinger",
    colour_var = "group",
    shape_var = "None",
    ellipse = FALSE
  )

  expect_s3_class(res$plot, "ggplot")
  expect_match(res$plot$labels$subtitle, "aitchison distance")
  expect_match(res$plot$labels$subtitle, "clr")
  expect_equal(nrow(res$data), phyloseq::nsamples(ps))
})
