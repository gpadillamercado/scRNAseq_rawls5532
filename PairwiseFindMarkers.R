# Pairwise Find Markers Loop
# This function allows pair-wise comparisons to be made and collated from the
# FindMarkers function in Seurat
# object is the Seurat object

# This returns a dataframe of all pairwise comparisons
PairwiseFindMarkers <- function(object, ...) {
  nclusters <- max(
    as.numeric(
      levels(
        object@meta.data$seurat_clusters)))
  datalist <- list()
  counter <- 1
  for (i in 0:nclusters) {
    j <- i + 1
    while (j <= nclusters) {
      print(paste0("Comparing clusters ",
                   levels(object@active.ident)[[i+1]],
                   " & ", levels(object@active.ident)[[j+1]]))
      temp <- tibble::rownames_to_column(
        Seurat::FindMarkers(object,
                            ident.1 = levels(object@active.ident)[[i+1]],
                            ident.2 = levels(object@active.ident)[[j+1]], ...),
        var = "gene")

      tempA <- dplyr::mutate(dplyr::filter(temp, avg_log2FC < 0),
                             avg_log2FC = -1 * avg_log2FC,
                             cluster1 = levels(object@active.ident)[[j+1]],
                             cluster2 = levels(object@active.ident)[[i+1]])
      datalist[[counter]] <- tempA
      counter <- counter + 1

      tempB <- dplyr::mutate(dplyr::filter(temp, avg_log2FC > 0),
                             avg_log2FC = avg_log2FC,
                             cluster1 = levels(object@active.ident)[[i+1]],
                             cluster2 = levels(object@active.ident)[[j+1]])
      datalist[[counter]] <- tempB
      counter <- counter + 1

      j <- j + 1
    }
  }
  AUdf <- dplyr::bind_rows(datalist)
  AUdf
}
