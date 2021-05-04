# Find All Differentially Expressed Genes between two conditions
#
# This function allows comparisons across conditions to be made and collated from the
# FindMarkers function in Seurat
# object is the Seurat object
# condition1 is the experimental condition
# condition2 is the second condition, usually a control

# Returns a combined dataframe of all comparisons

ConditionDiffExp <- function(object, condition1, condition2 = NULL, ...) {
  datalist <- list() # initialize an empty list where results will be stored
  nclusters <- max(
    as.numeric(
      levels(
        object@meta.data$seurat_clusters))) # store number of clusters in 'object'
  i <- 0
  while(i <= nclusters){ # loop ends when all clusters have had condition1 vs condition2 comparisons
    print(paste0("Comparing ", condition1," / ",condition2, "in cluster ",
                 levels(object@active.ident)[[i+1]])) # This prints the named identity of the cluster
    
    temp <- dplyr::mutate(tibble::rownames_to_column( # makes a temporary data.frame that will be stored in datalist
      Seurat::FindMarkers(object = object,
                          ident.1 = condition1,
                          ident.2 = condition2,
                          group.by = "library",
                          subset.ident = levels(object@active.ident)[[i+1]]), # logfc.threshold = 0, min.pct = 0.01
      var = "gene"),
      cluster = levels(object@active.ident)[[i+1]])
    datalist[[i+1]] <- temp # store results in datalist
    i <- i+1 
  }
  DEC <- dplyr::bind_rows(datalist) # combines all results from datalist into a single data.frame
  DEC
}
