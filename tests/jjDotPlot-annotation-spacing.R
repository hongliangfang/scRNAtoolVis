# A top marker-group annotation should be close to the panel even with
# many Seurat clusters, and an explicit y position must remain supported.

n_clusters <- 18L
counts <- matrix(rep(c(3, 2, 1), n_clusters), nrow = 3,
                 dimnames = list(c("GeneA", "GeneB", "GeneC"),
                                 paste0("cell", seq_len(n_clusters))))
object <- Seurat::CreateSeuratObject(
  counts = counts, min.cells = 0, min.features = 0
)
object$seurat_clusters <- as.character(seq_len(n_clusters) - 1L)
object <- Seurat::NormalizeData(object, verbose = FALSE)

marker_df <- data.frame(
  gene = c("GeneA", "GeneB", "GeneC"),
  cluster = c("TypeA", "TypeB", "TypeB")
)
plot_args <- list(
  object = object, assay = "RNA", slot = "data",
  id = "seurat_clusters", markerGene = marker_df,
  anno = TRUE, scale = FALSE, ytree = FALSE
)

# Each annotation_custom layer stores its y coordinate in geom_params.
annotation_ys <- function(p) {
  unlist(lapply(p$layers, function(layer) layer$geom_params$ymin))
}

p <- do.call(scRNAtoolVis::jjDotPlot, plot_args)
stopifnot(
  inherits(p, "ggplot"),
  any(abs(annotation_ys(p) - (n_clusters + 0.9)) < 1e-8)
)

custom_y <- n_clusters + 2
p_custom <- do.call(scRNAtoolVis::jjDotPlot, c(
  plot_args, list(anno_y_position = custom_y)
))
stopifnot(any(abs(annotation_ys(p_custom) - custom_y) < 1e-8))

# Numeric IDs use continuous ggplot y scaling rather than discrete levels.
object$seurat_clusters <- seq_len(n_clusters) - 1L
plot_args$object <- object
p_numeric <- do.call(scRNAtoolVis::jjDotPlot, plot_args)
expected_numeric_y <- (n_clusters - 1) * 1.05 + 0.3
stopifnot(any(abs(annotation_ys(p_numeric) - expected_numeric_y) < 1e-8))
