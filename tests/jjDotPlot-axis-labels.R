# Regression checks for x-axis labels when jjDotPlot rotates gene names.
# This script is run by R CMD check after the package has been installed.

counts <- matrix(
  c(5, 0, 3, 1,
    2, 1, 0, 6,
    0, 4, 5, 0),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("cell", 1:4))
)

object <- Seurat::CreateSeuratObject(
  counts = counts, min.cells = 0, min.features = 0
)
object$seurat_clusters <- factor(c("0", "0", "1", "1"))
object <- Seurat::NormalizeData(object, verbose = FALSE)

plot_args <- list(
  object = object, assay = "RNA", slot = "data",
  id = "seurat_clusters", gene = rownames(counts),
  scale = FALSE, ytree = FALSE
)

# A 45-degree rotation must bring gene labels close to their tick marks.
p45 <- do.call(scRNAtoolVis::jjDotPlot, c(plot_args, list(x.text.angle = 45)))
stopifnot(
  inherits(p45, "ggplot"),
  isTRUE(all.equal(p45$theme$axis.text.x$angle, 45)),
  isTRUE(all.equal(p45$theme$axis.text.x$hjust, 1)),
  isTRUE(all.equal(p45$theme$axis.text.x$vjust, 1))
)

# The default vertical layout must stay backward compatible.
p90 <- do.call(scRNAtoolVis::jjDotPlot, plot_args)
stopifnot(
  isTRUE(all.equal(p90$theme$axis.text.x$angle, 90)),
  isTRUE(all.equal(p90$theme$axis.text.x$vjust, 0))
)

# The caller must still be able to override automatic alignment.
p_custom <- do.call(scRNAtoolVis::jjDotPlot, c(
  plot_args, list(x.text.angle = 45, x.text.vjust = 0.5)
))
stopifnot(isTRUE(all.equal(p_custom$theme$axis.text.x$vjust, 0.5)))
