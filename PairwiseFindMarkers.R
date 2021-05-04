# Pairwise Find Markers Loop
# This function allows pair-wise comparisons to be made and collated from the
# FindMarkers function in Seurat
# object is the Seurat object

# This returns a dataframe of all pairwise comparisons
PairwiseFindMarkers <- function(object, ...) {
  nclusters <- max(
    as.numeric(
      levels(
        object@meta.data$seurat_clusters))) # stores number of clusters in 'object'
  datalist <- list() # initializes empty list for storing results
  counter <- 1 # initializes a counter to forgo calculating an index based on nested loops
  for (i in 0:nclusters) { # loop running for each clusters
    j <- i + 1
    while (j <= nclusters) { # loop running for each cluster that is "after" cluster 'i'
      print(paste0("Comparing clusters ",
                   levels(object@active.ident)[[i+1]],
                   " & ", levels(object@active.ident)[[j+1]])) # prints which clusters are being compared
      temp <- tibble::rownames_to_column( # create data.frame to store results
        Seurat::FindMarkers(object,
                            ident.1 = levels(object@active.ident)[[i+1]],
                            ident.2 = levels(object@active.ident)[[j+1]], ...),
        var = "gene")

      tempA <- dplyr::mutate(dplyr::filter(temp, avg_log2FC < 0), # is there decreased relative expression?
                             avg_log2FC = -1 * avg_log2FC,        # then reverse the signs of avg_log2FC
                             cluster1 = levels(object@active.ident)[[j+1]], # and switch the cluster comparison order
                             cluster2 = levels(object@active.ident)[[i+1]])
      datalist[[counter]] <- tempA # store this results in datalist
      counter <- counter + 1 # increase counter to next index in datalist

      tempB <- dplyr::mutate(dplyr::filter(temp, avg_log2FC > 0), # if there is increased expression
                             avg_log2FC = avg_log2FC,             # keep the sign of avg_log2FC and cluster comparison order
                             cluster1 = levels(object@active.ident)[[i+1]],
                             cluster2 = levels(object@active.ident)[[j+1]])
      datalist[[counter]] <- tempB # store this results in datalist
      counter <- counter + 1 # increase counter to next index in datalist

      j <- j + 1
    } # once this loop runs for all clusters 'after' the current cluster, it bumps up to the 'parent' loop
  }   # the changing of cluster order and sign of log2FC reduces FindMarkers computations from (n * (n-1)) to (n * (n-1) - (n * (n-1))/2)
      # but is only looking at increased expression (which makes sense if you are looking for markers of clusters)
  
  AUdf <- dplyr::bind_rows(datalist) # combine all results into a single data.frame
  AUdf
}
