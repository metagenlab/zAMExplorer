test_that("publication theme and palettes are reusable", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) + ggplot2::geom_point() + zAMPExplorer:::zamp_theme()
  expect_s3_class(p, "ggplot"); expect_length(zAMPExplorer:::zamp_palette(12), 12)
})

test_that("plot export writes a PDF", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) + ggplot2::geom_point(); file <- tempfile(fileext = ".pdf")
  zAMPExplorer:::zamp_export_plot(p, file, "pdf", 4, 3); expect_true(file.exists(file)); expect_gt(file.info(file)$size, 0)
})
