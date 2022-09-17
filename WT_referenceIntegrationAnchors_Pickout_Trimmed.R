#Package and Data Loading----

library(Seurat)
library(tidyverse)
library(DESeq2)

# Read the file in for quick analyses
df_integrated <- readRDS("Raw_Data/rawls5532_reference_integrated_GPM_UMAP52_35DimIntRes0.82.rds")
df_integrated.cluster17 <- readRDS(file = "scData/rawls5532_refFinal_GPM_17_2clusters_res0.5.rds")

DefaultAssay(df_integrated.cluster17) <- "RNA"

# Integrating based on a reference dataset
df_2 <- readRDS("scData/rawls_5532_Seuratv3_1.Robj")

WT <- readRDS("scData/WT_original_gpm.rds")
MUT <- readRDS("scData/MUT_original_gpm.rds")

# Integration of Datasets----
# Create a variable which specifies genotype
WT$library <- "WT"
MUT$library <- "MUT" 

# Make a list of these object
df_list <- list("WT" = WT,"MUT" = MUT)
for (i in 1:length(df_list)) {
  # Create a feature that finds proportion of mitochondrial genes in a cell
  df_list[[i]] <- PercentageFeatureSet(df_list[[i]],
                                               pattern = "^mt-", 
                                               col.name = "percent.mt")
}
for (i in 1:length(df_list)) {
  df_list[[i]] <- subset(df_list[[i]],
                                 subset = nCount_RNA < 30000 & 
                                   nFeature_RNA < 3000 & 
                                   percent.mt < 25)
}

# Performs Normalization on the list of Seurat objects
df_list <- lapply(X = df_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
PC.int <- 35
# Selects the features to integrate by:
features <- SelectIntegrationFeatures(object.list = df_list, dims = 1:PC.int)

# Runs PCA for rpca reduction
df_list <- lapply(X = df_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# Get the WT dataset.
reference_dataset <- which(names(df_list) == "WT")

# Find the Anchors using the reference
df_anchors <- FindIntegrationAnchors(object.list = df_list,
                                     reference = reference_dataset,
                                     reduction = "rpca",
                                     dims = 1:PC.int,
                                     verbose = TRUE)

df_integrated <- IntegrateData(anchorset = df_anchors, dims = 1:PC.int)

# PCA and Jackstraw----

# Do PCA on this integrated data
DefaultAssay(df_integrated) <- "integrated"
df_integrated <- ScaleData(object = df_integrated, verbose = FALSE)
df_integrated <- RunPCA(object = df_integrated,
                                features = VariableFeatures(object = df_integrated),
                                npcs = 100)

ElbowPlot(df_integrated, ndims = 100)
# Redo Jackstraw for the reference-integrated dataset
df_integrated <- JackStraw(object = df_integrated,
                           num.replicate = 100,
                           dims = 100,
                           prop.freq = 0.1,
                           verbose = TRUE)
df_integrated <- ScoreJackStraw(object = df_integrated,
                                dims = 1:100) 
#
JackStrawPlot(object = df_integrated, dims = 1:60)

# Chose 52 PCs to Run the UMAP algorithm based on elbow and Jackstraw Plot
# UMAP----
pcs.use <- 52 # The selected number of clusters
df_integrated2 <- RunUMAP(df_integrated, 
                           reduction = "pca", dims = 1:pcs.use)
df_integrated2 <- FindNeighbors(df_integrated2, 
                                 reduction = "pca", dims = 1:pcs.use)


df_integrated2 <- FindClusters(df_integrated2, 
                              resolution = 0.82)

DimPlot(object = df_integrated2, reduction = "umap", split.by = "library",
        label = TRUE)
ggsave("WTReferenced52PCs_35DimInts_Res0.82.pdf", height = 5, width = 11)

# for (i in 0:4) {
#   df_integrated <- FindClusters(df_integrated, 
#                                 resolution = 0.6+(i/20))
#   temp <- DimPlot(object = df_integrated, reduction = "umap",
#         split.by = "library", label = TRUE)
#   png(filename=paste0("scData/OUTPUT/MUT_referencedUMAP_52DimsRes",
#                       0.6+(i/20),
#                       "_50intDims.png"),
#       res=400,width=4800,height=2400)
#   print(temp)
#   dev.off()
#   print(i)
# }

##
# Painting----
# Define the number of clusters in the source clustering map
SourceCluster <- 22

# Makes an empty list to put values
vectorClusters <- list()
# Reads all the UMIs from MUT.UMAP
allUMI.UMAP <- colnames(df_2)

# This subsets to only those UMIs in cluster 8 of the integrated analysis
# The UMI with a '_1' attached
# All the UMIs from the WT only analysis (no "_1" to extract)
# Returns the indexes in allUMI_MUT.UMAP that match for each selected UMI list
for (i in 0:SourceCluster) {
  print(i)
  ClusterVector <- WhichCells(object = df_2, idents = i) %>%
    base::match(allUMI.UMAP)
  vectorClusters[[i+1]] <- ClusterVector
}

# Put these values in a dataframe for easy access/organization
vectorClusters2 <- tibble(Idents = c(0:SourceCluster)) %>%
  mutate(DefClust = vectorClusters)

UMAP_painted <- df_integrated2

for (i in 0:SourceCluster){
  Idents(UMAP_painted,
         cells = vectorClusters2$DefClust[[i+1]]) <- vectorClusters2$Idents[i+1]
}
levels(UMAP_painted) <- c(0:22)
table(Idents(UMAP_painted))
table(Idents(df_2))


# Assign ggplots as objects for cowplot output
tLp <- DimPlot(UMAP_painted, label = TRUE) + ggtitle("Transferred Labels")
oLp <- DimPlot(df_integrated2, label = TRUE) + ggtitle("New Plot, FINAL Clustering")
rp <- DimPlot(df_2, label = TRUE) + ggtitle("Original Plot, from Core")

cowplot::plot_grid(rp, tLp, oLp, nrow = 1)
ggsave("OriginalClusterLabels_Transferred2FinalUMAP.pdf", height = 5, width = 17)
##

DefaultAssay(df_integrated2) <- "RNA"
# Normalize RNA data for visualization purposes
df_integrated2 <- NormalizeData(df_integrated, verbose = FALSE)
# Save the object at this point 
saveRDS(df_integrated2, "scData/rawls5532_reference_integrated_GPM_UMAP52_35DimIntRes0.82.rds")

# Read the file in
df_integrated <- readRDS("scData/rawls5532_reference_integrated_GPM_UMAP52_35DimIntRes0.82.rds")

DimFeaturePlot(df_integrated2,
            features = c('fabp2', 'nr1h4','slc10a2'), split.by = 'library', pt.size = 0.1)
ggsave("scData/WT_referenced60_Normalized_GeneSet1.pdf", height = 8, width = 4)

# Looking at all fabp6 expressing cells
PickingOut_fabp6 <- WhichCells(object = df_integrated,
                               expression = fabp6 > 2) %>%
  base::match(colnames(df_integrated))
rightPlot <- DimPlot(df_integrated, cells = PickingOut_fabp6,
                     group.by = "library", pt.size = 0.8)
leftPlot <- DimPlot(df_integrated, cells = PickingOut_fabp6,
                    label = TRUE, pt.size = 0.8)
cowplot::plot_grid(leftPlot, rightPlot)
ggsave("Pick_fabp6.pdf", height = 5, width = 11)


# Find Conserved Markers
conMarkers <- FindConservedMarkers(df_integrated, ident.1 = 10,
                                   grouping.var = "library")

VlnPlot(df_integrated, features = c("fabp6","slc10a2"), split.by = "library",
        assay = "RNA")

df_integrated3 <- subset(df_integrated, subset = orig.ident == "MUT")
# Find All the Markers for each cluster
refMarkers <- FindAllMarkers(object = df_integrated,
                             only.pos = TRUE,
                             min.pct = 0.1,
                             thresh.use = 0.25,
                             assay = "RNA")
top50 <- refMarkers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

write_csv(top50, "scData/WTReferencedWTonly_top50_Markers_52PCs35DimIntsres0.8.csv")
write_csv(refMarkers, "scData/WTReferenced_Markers_52PCs35DimIntsres0.8.csv")

###### Generating Figures for the Paper Here ----
# These are draft figures
markers <- read_csv("scData/WTReferencedCombined_top50_Markers_52PCs35DimIntsres0.8.csv")
markers5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
markers100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DotPlot(df_integrated, features = unique(markers5$gene), dot.scale = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("PaperFigures/Figure2_MarkerGenesDotPlot_20200305.pdf", height = 6, width = 18)

DimPlot(subset(df_integrated, subset = orig.ident == "WT"), label = TRUE)
ggsave("PaperFigures/Figure1_WT.UMAP_20200305.pdf", height = 8, width = 10)

FeaturePlot(subset(df_integrated, subset = orig.ident == "WT"), 
            features = c("fabp2", "agr2", "slc10a2", "fabp6", "her6", "neurod1"))
ggsave("PaperFigures/Figure2_WT.UMAP_features_20200305.pdf", height = 14, width = 10)


DoHeatmap(df_integrated, features = markers100$gene, assay = "integrated") + 
  NoLegend() +
  theme(axis.text.y = element_text(size = 8))
ggsave("PaperFigures/FigureS1_featuresHeatmap_20200305.pdf", height = 14, width = 20)
####

# Quick MUT v WT FindMarkers for ALL clusters----
for (i in 0:26L) {
  FindMarkers(object = df_integrated,
              ident.1 = "MUT",
              group.by = "library",
              subset.ident = i) %>%
    rownames_to_column(var = "gene") %>% 
    top_n(n = 500, wt = avg_logFC) %>%
    write_csv(paste0("cluster", i, "_diff", ".csv"))
}

FindMarkers(object = df_integrated, ident.1 = c(0, 2, 15)) %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("clusters", "0_2_15", "_markers", ".csv"))

FindMarkers(object = df_integrated, ident.1 = c(1, 7, 10)) %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("clusters", "1_7_10", "_markers", ".csv"))

FindMarkers(object = df_integrated, ident.1 = 5) %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("clusters", "5", "_markers", ".csv"))

FindMarkers(object = df_integrated, ident.1 = 6) %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("clusters", "6", "_markers", ".csv"))

FindMarkers(object = df_integrated, ident.1 = c(8, 11, 12, 21, 22)) %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("clusters", "8_11_12_21_22", "_markers", ".csv"))

FindMarkers(object = df_integrated, ident.1 = c(0, 15, 2), ident.2 = c(1, 10, 7)) %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("clusters", "0_15_2_vs_1_7_10", "_markers", ".csv"))

# Exploratory data analysis of mut v wt----
library(readxl)    
read_excel_allsheets <- function(filename, tibble = TRUE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X) %>% 
                filter(gene %in% markers$gene))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

MvW <- read_excel_allsheets("01022020_scRNASeq_Mut_vs_WT_EachCluster_top500.xlsx")

DotPlot(df_integrated2,
        features = top5$gene %>% unique(),
        split.by = 'library', cols = c('blue','red')) + RotatedAxis()


FeaturePlot(df_integrated2, features = "percent.mt", split.by = "library")

DotPlot(df_integrated2,
        features = c('fabp6', 'slc10a2',
                     'dab2', 'cubn', 'amn', 'tmigd1'),
        split.by = "library", cols = c("blue", "red"))

DimPlot(df_integrated2, label = TRUE)
FeaturePlot(df_integrated2, features = c('elf3'))

df_int2Avg <- AverageExpression(df_integrated2, assays = "RNA")
write.csv(df_int2Avg, "AverageGeneExpression_FinalCluster.csv")

sdf <- subset(df_integrated, subset = orig.ident == "WT")
df_int2Avg <- AverageExpression(sdf, assays = "RNA")
write.csv(df_int2Avg, "AverageGeneExpression_Wild-type_FinalCluster.csv")

sdf_m <- subset(df_integrated,subset = orig.ident == "MUT")
df_int2Avg <- AverageExpression(sdf_m, assays = "RNA")
write.csv(df_int2Avg, "AverageGeneExpression_fxr-Mutant_FinalCluster.csv")

library(clipr)
write_clip(table(Idents(df_integrated)))
write_clip(table(Idents(sdf_m)))

##
# Do a pseudobulk analysis----
##
WTpseudobulk <- as.data.frame(WT@assays$RNA@counts) %>% rownames_to_column(var = "gene")
MUTpseudobulk <- as.data.frame(MUT@assays$RNA@counts) %>% rownames_to_column(var = "gene")

length(colnames(WTpseudobulk)) # There are 5661 WT cells
length(colnames(MUTpseudobulk)) # There are 6881 MUT cells


MUTpseudoBulk <- MUTpseudobulk %>% filter(gene %in% WTpseudobulk$gene) %>%
  column_to_rownames(var = 'gene')
WTpseudoBulk <- WTpseudobulk %>% filter(gene %in% MUTpseudobulk$gene) %>%
  column_to_rownames(var = 'gene')
# Logic of bulk gene expression vector

log10(1+((10^4)*(apply(EEdf, 1, sum)/sum(EEdf))))
# So the gene expression vector is the sum of all the counts for each gene
# Then I will divide this by the total number UMIs for this sample (WT or MUT)
# Then multiply by 10^4 and add 1
# Finally take the log10 of this value

WTpseudobulkRes <- log10(1+((10^4)*(apply(WTpseudobulk, 1, sum)/sum(WTpseudobulk))))
MUTpseudobulkRes <- log10(1+((10^4)*(apply(MUTpseudobulk, 1, sum)/sum(MUTpseudobulk))))

CombopseudoBulk1 <- tibble(Gene =  names(WTpseudobulkRes),
                          WT =  WTpseudobulkRes)
CombopseudoBulk2 <- tibble(Gene =  names(MUTpseudobulkRes),
                           MUT =  MUTpseudobulkRes)
CombopseudoBulk <- inner_join(CombopseudoBulk1, CombopseudoBulk2)


ggplot(CombopseudoBulk, aes(x = WT, y = MUT)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") + theme_bw() +
  ggtitle("Condition-Normalized Pseudobulk RNA-Seq")
ggsave("Pseudobulk_combination.pdf", height = 5, width = 5)

longCombopseudoBulk <- CombopseudoBulk %>% 
  pivot_longer(2:3, names_to = "Sample", values_to = "Expression")
levels(longCombopseudoBulk) <- c('WT', 'MUT')

geneset <- c("fabp2", "fabp6", "neurod1", "nr1h4", "slc10a2")

ggplot(data = (longCombopseudoBulk %>% filter(Gene %in% geneset)),
       aes(x = Gene, y = Expression, group = Gene, color = Sample)) +
  geom_point(position = position_dodge(width = 0.75)) + theme_bw() +
  ggtitle("Average Normalized Expression from Pseudobulk scRNA-Seq data")


##
# Subsample and putting into DESeq2
Sample_1 <- floor(runif(length(colnames(WTpseudoBulk))*0.6,
      min = 1, max = 1+length(colnames(WTpseudoBulk))))  
Sample_2 <-  floor(runif(length(colnames(WTpseudoBulk))*0.6,
                 min = 1, max = 1+length(colnames(WTpseudoBulk))))  
Sample_3 <-  floor(runif(length(colnames(WTpseudoBulk))*0.6,
                 min = 1, max = 1+length(colnames(WTpseudoBulk))))  

# Subsample and putting into DESeq2
Sample_1m <- floor(runif(length(colnames(MUTpseudoBulk))*0.6,
                        min = 1, max = 1+length(colnames(MUTpseudoBulk))))  
Sample_2m <-  floor(runif(length(colnames(MUTpseudoBulk))*0.6,
                         min = 1, max = 1+length(colnames(MUTpseudoBulk))))  
Sample_3m <-  floor(runif(length(colnames(MUTpseudoBulk))*0.6,
                         min = 1, max = 1+length(colnames(MUTpseudoBulk))))  

cts <- data.frame(
  gene_id = row.names(MUTpseudoBulk),
  MUT_1 = MUTpseudoBulk %>% dplyr::select(Sample_1m) %>% base::apply(1, sum),
  MUT_2 = MUTpseudoBulk %>% dplyr::select(Sample_2m) %>% base::apply(1, sum),
  MUT_3 = MUTpseudoBulk %>% dplyr::select(Sample_3m) %>% base::apply(1, sum),
  WT_1 =  WTpseudoBulk %>% dplyr::select(Sample_1) %>% base::apply(1, sum),
  WT_2 =  WTpseudoBulk %>% dplyr::select(Sample_2) %>% base::apply(1, sum),
  WT_3 =  WTpseudoBulk %>% dplyr::select(Sample_3) %>% base::apply(1, sum)
)

coldata <- data.frame(
  c("MUT_1","MUT_2","MUT_3",
    "WT_1","WT_2","WT_3"),
  library = c("MUT","MUT","MUT",
              "WT","WT","WT"),
  batch = c("1", "2", "3","1", "2", "3"))
write_tsv(coldata, "colData.txt")
coldata <- read.csv("colData.txt",
                sep = "\t",
                row.names = 1)
head(coldata)

# DESeq2 with pseudobulked values----
# Most of this was redundant with previous scRNA analysis but I did not realize at the time

library(DESeq2)

write_tsv(cts, "SampledSCRNA_full.txt")
cts <- as.matrix(read.csv("SampledSCRNA_full.txt",
                          sep = "\t",
                          row.names = "gene_id"))
head(cts)

# Load in the sampled count totals
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ library
)

dds$library <- relevel(dds$library, "WT")
ddsWald <- DESeq(dds)
library(ashr)


resultsNames(ddsWald)
                          

subResults <- results(ddsWald, name = "library_MUT_vs_WT")

plot(subResults$log2FoldChange, -1*log10(subResults$padj))

ddsWald_shrunk <- lfcShrink(ddsWald, contrast = c("library", "MUT", "WT"),
                            res = subResults, type = "ashr")

plot(ddsWald$log2FoldChange, -1*log10(ddsWald$padj),
     xlab =expression(log[2]*"(MUT/WT), FXR"),
     ylab =expression(-log[10]*"(Adjusted p value)"),
     xlim =c(-10, 10),
     main = "Volcano Plot Pseudobulk Subsamples N = 3 per Genotype",
     pch = 16,
     cex = 0.5,
     col =rgb(0, 0, 0))

head(ddsWald_shrunk)


resdf <- as.data.frame(ddsWald_shrunk) %>% rownames_to_column("gene")
resdf_filtered <- as.data.frame(ddsWald_shrunk) %>% rownames_to_column("gene") %>%
  filter(padj <= 0.01 & (log2FoldChange >= 1 | log2FoldChange <= -1))

write_csv(resdf, "Pseudobulked_scRNAseq_DESeq2_60percentWT_3reps_total.csv")
write_csv(resdf_filtered, "FilteredPseudobulked_scRNAseq_DESeq2_60percentWT_3reps.csv")
write_csv(CombopseudoBulk %>% filter(((MUT/WT) >= 2 | (MUT/WT) <= 0.5) & is.finite(MUT/WT)),
          "Peudobulked_scRNAseq_Normalized.csv")

##
## Misc plots
##
Mut <- readRDS("MUT_UMAP.rds")
DimPlot(Mut, label = TRUE) + ggtitle("Mutant-only clustered with 40 PCs")
ggsave("scData/MutantUMAP_labelled.pdf", height = 6, width = 7.5)

Wt <- readRDS("scData/WT_UMAP.rds")
DimPlot(Wt, label = TRUE) + ggtitle("Wildtype-only clustered with 37 PCs")
ggsave("scData/WTtUMAP_labelled.pdf", height = 6, width = 7.5)


DimPlot(df, label = TRUE) + ggtitle("Integrated data (from core) clustered with 52 PCs")
ggsave("scData/IntegratedUMAP_labelled.pdf", height = 6, width = 7.5)

df_integrated <- readRDS("scData/rawls5532_reference_integrated_GPM_UMAP.rds")
bottom <- DimPlot(df_integrated, label = TRUE, split.by = "library") +
  ggtitle("WT-referenced integrated data clustered with 52 PCs, by genotype")
ggsave("scData/WTrefd_IntegratedUMAP_labelled.pdf", height = 6, width = 15)
top <- DimPlot(df, label = TRUE, split.by = "library") +
  ggtitle("Integrated data clustered with 52 PCs, by genotype")
ggsave("scData/WTrefd_IntegratedUMAP_labelled_split.pdf", height = 6, width = 15)


plot_grid(top, bottom, nrow = 2)
ggsave("scData/Both_IntegratedUMAP_labelled_split.pdf", height = 12, width = 15)

df.no.umap <- df_integrated
df.no.umap[["umap"]] <- NULL
FeaturePlot(df.no.umap, feature = 'elf3', label = TRUE) + RotatedAxis()
FeaturePlot(df_integrated, feature = 'elf3', label = TRUE) + RotatedAxis()

FeatureScatter(df_integrated, 'fabp6', 'slc10a2', group.by = "library")

##
## Subsetting cluster 17----
##
df_integrated.cluster17 <- subset(df_integrated, idents = '17')
DefaultAssay(df_integrated.cluster17) <- "integrated"
df_integrated.cluster17 <- FindVariableFeatures(df_integrated.cluster17,
                                                verbose = FALSE)
df_integrated.cluster17 <- ScaleData(df_integrated.cluster17, verbose = FALSE)
df_integrated.cluster17 <- RunPCA(object = df_integrated.cluster17,  npcs = 100)
df_integrated.cluster17 <- ProjectDim(object = df_integrated.cluster17)

df_integrated.cluster17 <- JackStraw(object = df_integrated.cluster17,
                                     num.replicate = 100,
                                     dims = 100,
                                     prop.freq = 0.1, verbose = TRUE)
df_integrated.cluster17 <- ScoreJackStraw(object = df_integrated.cluster17,
                                          dims = 1:100) 
JackStrawPlot(object = df_integrated.cluster17, dims = 1:100)

pcs.use=8

df_integrated.cluster17 <- FindNeighbors(object = df_integrated.cluster17,
                                         dims = 1:pcs.use)
df_integrated.cluster17 <- FindClusters(object = df_integrated.cluster17,
                                        resolution = 0.5)
df_integrated.cluster17 <- RunUMAP(object = df_integrated.cluster17,
                                   dims = 1:pcs.use,
                                   assay = "integrated")
DimPlot(df_integrated.cluster17) + ggtitle(label = "Sub-clustering Cluster 17")

subclust_0 <- WhichCells(object = df_integrated.cluster17,
                         idents = "0") %>%
  base::match(colnames(df_integrated))
subclust_1 <- WhichCells(object = df_integrated.cluster17,
                     idents = "1") %>%
  base::match(colnames(df_integrated))
subclust_2 <- WhichCells(object = df_integrated.cluster17,
                         idents = "2") %>%
  base::match(colnames(df_integrated))

# Looking at all fabp6 expressing cells

leftPlot <- DimPlot(df_integrated, cells.highlight = subclust_0,
                     group.by = "library",
                     label = TRUE, pt.size = 0.8) + ggtitle("Sub-cluster 0") +
  theme(legend.position = "none")
midPlot <- DimPlot(df_integrated, cells.highlight = subclust_1, 
                    group.by = "library",
                    label = TRUE, pt.size = 0.8) + ggtitle("Sub-cluster 1") +
  theme(legend.position = "none")
rightPlot <- DimPlot(df_integrated, cells.highlight = subclust_2, 
                     group.by = "library",
                     label = TRUE, pt.size = 0.8) + ggtitle("Sub-cluster 2") +
  theme(legend.position = "none")
cowplot::plot_grid(leftPlot, midPlot, rightPlot, align = "h")
ggsave("Cluster17_3SubclustersInContext.pdf", height = 5, width = 11)

# Set Identifiers
Idents(object = df_integrated)
Idents(object = df_integrated, cells = subclust_0) <- '17_0'
Idents(object = df_integrated, cells = subclust_1) <- '17_1'
Idents(object = df_integrated, cells = subclust_2) <- '17_2'
table(Idents(object = df_integrated))

DimPlot(df_integrated, label = TRUE, split.by = "library")

saveRDS(object = df_integrated, file = "scData/rawls5532_refFinal_GPM_17subclustered_res0.5.rds")
df_integrated_sub <- readRDS(file = "scData/rawls5532_refFinal_GPM_17subclustered_res0.5.rds")
saveRDS(object = df_integrated.cluster17, file = "scData/rawls5532_refFinal_GPM_17_2clusters_res0.5.rds")
df_integrated.cluster17 <- readRDS(file = "scData/rawls5532_refFinal_GPM_17_2clusters_res0.5.rds")

# Finding the differentiating markers
FindMarkers(object = df_integrated, ident.1 = "17_0") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_2subclusterRes0.5_", "0", "_markers", ".csv")) # Finding the differentiating markers
FindMarkers(object = df_integrated, ident.1 = "17_1") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_2subclusterRes0.5_", "1", "_markers", ".csv"))
# Finding the differentiating markers
FindMarkers(object = df_integrated, ident.1 = "17_2") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_3subcluster", "2", "_markers", ".csv"))
# Finding the differentiating markers between just 17 subclusters
FindMarkers(object = df_integrated, ident.1 = "17_1", ident.2 = "17_0") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_2subclusterRes0.5_", "1vs0", "_markers", ".csv"))

FindMarkers(object = df_integrated,
            ident.1 = "MUT",
            group.by = "library",
            subset.ident = "17_1") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_2subclusterRes0.5_", "1", "_MUTvWTdiff", ".csv"))

FindMarkers(object = df_integrated,
            ident.1 = "MUT",
            group.by = "library",
            subset.ident = "17_0") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_2subclusterRes0.5_", "0", "_MUTvWTdiff", ".csv"))

FindMarkers(object = df_integrated,
            ident.1 = "MUT",
            group.by = "library",
            subset.ident = "17_2") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("Clust17_3subcluster", "2", "_MUTvWTdiff", ".csv"))

DefaultAssay(df_integrated.cluster17) <- "RNA"

FeaturePlot(object = df_integrated,
            features = c("scpep1", "lamp2"),
            split.by = "library")
FeaturePlot(object = df_integrated.cluster17,
            features = c("fabp6", "amn"),
            split.by = "library",
            label = TRUE)
DimPlot(df_integrated, label = TRUE)
ggsave("PaperFigures/Figure4_FullUMAP_Cluster17_2subclusters_20200305.pdf", height = 8, width = 10)
clust17top10 <- FindMarkers(object = df_integrated.cluster17, ident.1 = 1) %>% rownames_to_column("gene") %>% write_csv("Cluster17_2subcluster_findmarkers_1.csv")
DotPlot(df_integrated.cluster17, features = clust17top10$gene,
        split.by = "library", cols = c("blue", "red")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("PaperFigures/Figure4_DotPlottop50_Cluster17_2subclusters_20200305.pdf", height = 4, width = 10)
table(Idents(subset(df_integrated.cluster17, subset = orig.ident == "MUT")))

### Differential Genes in each cluster
# Finding the differentiating markers
FindMarkers(object = df_integrated, ident.1 = "17_0") %>%
  rownames_to_column(var = "gene") %>% 
  write_csv(paste0("FinalClustering_", "0", "_markers", ".csv"))

df_2 <- df_2 %>% rownames_to_column(var = "Genes")

write_tsv(df_2, "STREAM_file.tsv")



table(Idents(df_integrated))

sigMarkers <- refMarkers %>% 
  filter(p_val_adj <= 0.05, !(avg_logFC = 0)) %>% group_by(cluster)

for (i in 0:26L) {
  FindMarkers(object = df_integrated,
              ident.1 = "MUT",
              group.by = "library",
              subset.ident = i,
              test.use = "DESeq2") %>%
    rownames_to_column(var = "gene") %>% mutate(Cluster = i) %>%
    write_csv(paste("MutvWT/DESeq2_Method/MutvWT_Cluster",
      i, "FinalUMAP_DESeq2.csv", sep = "_"))
}
# Then I did `cat MutvWT_Cluster* > list.csv` in terminal inside
# This generated list.csv with a fully appended clusters
mvw_list <- read_csv("MutvWT/DESeq2_Method/list.csv")
# Went to Excel and sorted `Cluster` and removed all the redundant text rows
df_mvw <- read_csv("MutvWT/DESeq2_Method/MutvWT_AllClusters_FinalUMAP.csv")

df_mvw <- read_csv("MutvWT/MutvWT_AllClusters_FinalUMAP.csv")

df_mvw2 <- df_mvw %>% group_by(Cluster) %>%
  filter(p_val < 0.05, !(avg_logFC = 0)) %>%
  summarise(Total = n(),
            Genes_Up = sum(avg_logFC > 0),
            Genes_Down = sum(avg_logFC < 0)) %>%
  pivot_longer(cols = c(Genes_Up, Genes_Down))

ggplot(df_mvw2) + 
  geom_col(aes(x = as.factor(Cluster),
               y = value,
               fill = name)) + theme_bw() +
  labs(x = "Cluster", y = "Number of Genes", fill = "") +
  scale_fill_manual(values = c("skyblue", "purple2")) +
  theme(axis.text = element_text(size = 8, hjust = 0.5, angle = 45))
ggsave("DifferentialGenes_Absolute_take2.pdf", height = 4, width = 5.5)

ggplot(df_mvw2) + 
  geom_col(aes(x = as.factor(Cluster),
               y = value/Total,
               fill = name)) + theme_bw() +
  labs(x = "Cluster", y = "Number of Genes", fill = "") +
  scale_fill_manual(values = c("skyblue", "purple2")) +
  theme(axis.text = element_text(size = 8, hjust = 0.5, angle = 45))
ggsave("DifferentialGenes_Relative_take2.pdf", height = 4, width = 5.5)

df_mvw3 <- df_mvw %>% group_by(Cluster) %>%
  filter(!(avg_logFC = 0)) %>%
  summarise(Total = n(),
            Genes_Up = sum(avg_logFC > 0.5),
            Genes_Down = sum(avg_logFC < -0.5),
            filteredTotal = sum(Genes_Up + Genes_Down)) %>%
  pivot_longer(cols = c(Genes_Up, Genes_Down))

ggplot(df_mvw3) + 
  geom_col(aes(x = as.factor(Cluster),
               y = value/filteredTotal,
               fill = name)) + theme_bw() +
  labs(x = "Cluster", y = "Number of Genes", fill = "") +
  scale_fill_manual(values = c("skyblue", "purple2")) +
  theme(axis.text = element_text(size = 8, hjust = 0.5, angle = 45))
ggsave("DifferentialGenes_Relative_0.5logFCthreshold.pdf", height = 4, width = 5.5)



# Integrating Monocle3 into pseudotime analysis-----
# loading
library(monocle3)
clust.keep <- as.character(c(0:17, 19, 21:24, 26)) # IEC + esohageal + phayngeal only
clust.keep <- as.character(c(0:26)) # All of the clusters
IEC_only <- subset(df_integrated, idents = clust.keep)

# Gene annotations
gene_annotation <- as.data.frame(rownames(IEC_only@reductions[["pca"]]@feature.loadings),
  row.names = rownames(IEC_only@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# Cell Metadata
cell_metadata <- as.data.frame(IEC_only@assays[["RNA"]]@counts@Dimnames[[2]],
  row.names = IEC_only@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# Counts Sparse Matrix
new_matrix <- IEC_only@assays[["RNA"]]@counts
new_matrix <- new_matrix[rownames(IEC_only@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- new_matrix


# Constructing basci Monocle CDS object
cds_from_seurat <- new_cell_data_set(expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation)

# Create partitions
partsDf <- data.frame(Clusters = IEC_only@meta.data$seurat_clusters, 
  UMI = IEC_only@assays$RNA@counts@Dimnames[[2]])

# new.idents_grant <- c("Ionocytes 1", "Esophageal 1", "Ionocytes 2", "Goblet cells 2",
#   "Enterocytes (segment 1)", "preEEC", "preSec", "Esophageal 2", 
#   "Enteroendocrine 1", "Enterocytes (segment 3)",
#   "Esophageal 3", "Enteroendocrine 4", "Enteroendocrine 5", "Goblet cells 3", 
#   "Goblet cells 1", "Ionocytes 3", "Enterocytes (segment 1 precursors)",
#   "Enterocytes (segment 2)", "Neutrophils/Macrophages", 
#   "ISCs?", "Exocrine Pancrease", "Enteroendocrine 3", "Enteroendocrine 2",
#   "Pharyngeal 1", "Pharyngeal 2", "RBCs", "Pharyngeal 3")

partsDf <- partsDf %>%
  mutate(supercluster = case_when(
    Clusters %in% c(0, 2, 4, 9, 15, 16, 17) ~ "Enterocytes",
    Clusters %in% c(3, 6, 13, 14) ~ "Secretory",
    Clusters %in% c(5, 8, 11, 12, 21, 22) ~ "EECs",
    Clusters %in% c(1, 7, 10, 18, 19, 20, 23, 24, 25, 26) ~ "Non-IECs")
    ) %>%
  mutate(partitions = case_when(
    supercluster == "Enterocytes" ~ 1,
    supercluster == "Secretory" ~ 2,
    supercluster == "EECs" ~ 3,
    supercluster == "Non-IECs" ~ 4))

# add in the partitions
recreate.partition <- c(partsDf$partitions)
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

# Assign the cluster information
list_cluster <-  IEC_only@meta.data[["seurat_clusters"]]
names(list_cluster) <- IEC_only@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster


# Placeholder for louvian parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- NA

# Assign UMAP coordinates
cds_from_seurat@reducedDims@listData[["UMAP"]] <- IEC_only@reductions[["umap"]]@cell.embeddings

# Assign feature loading for downstream analysis
cds_from_seurat@preprocess_aux$gene_loadings <- IEC_only@reductions[["pca"]]@feature.loadings

# Learn graph (takes some time)
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions == "2"] <- "1"

head(partitions(cds_from_seurat))

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

# Plot this
plot_cells(cds_from_seurat,
           color_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 3)

plot_cells(cds_from_seurat,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
# Choose the root node?
# Someworkaround to color_cells_by "pseudotime"
## https://github.com/cole-trapnell-lab/monocle3/issues/130
## Do this before order_cells()
rownames(cds_from_seurat@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds_from_seurat@reducedDims$UMAP) <- NULL

# Set the root with GUI... this worked!!!
cds_from_seurat <- order_cells(cds_from_seurat)


# Plot by pseudotime
plot_cells(cds_from_seurat,
          group_cells_by = "partition",
          color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
  show_trajectory_graph = FALSE)

cds_subset <- choose_cells(cds_from_seurat)
subset_pr_test_res <- graph_test(cds_from_seurat, neighbor_graph = "principal_graph", cores = 4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_from_seurat[pr_deg_ids,], resolution=0.001)

cell_group_df <- tibble::tibble(cell = row.names(colData(cds_from_seurat)), 
                                cell_group = colData(cds_from_seurat))

agg_mat <- aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))

module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(cds_from_seurat,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)

someGenes <- c('gcga', 'insl5a', 'penka5a', 'trpa1b', 'pyyb')
someGenes_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% someGenes,]

plot_genes_in_pseudotime(someGenes_cds, color_cells_by = "clusters", min_expr = 0.5)

# -----------------------------------------
# Comparison with Bagnat/Sadler data
# Read the file in
df_integrated <- readRDS(
  "scData/rawls5532_reference_integrated_GPM_UMAP52_35DimIntRes0.82.rds")
new.idents_grant <- c("Ionocytes 1", "Esophageal 1", "Ionocytes 2", "Goblet cells 2",
  "Enterocytes (segment 1)", "preEEC", "preSec", "Esophageal 2", 
  "Enteroendocrine 1", "Enterocytes (segment 3)",
  "Esophageal 3", "Enteroendocrine 4", "Enteroendocrine 5", "Goblet cells 3", 
  "Goblet cells 1", "Ionocytes 3", "Enterocytes (segment 1 precursors)",
  "Enterocytes (segment 2)", "Neutrophils/Macrophages", 
  "ISCs?", "Exocrine Pancrease", "Enteroendocrine 3", "Enteroendocrine 2",
  "Pharyngeal 1", "Pharyngeal 2", "RBCs", "Pharyngeal 3")

new.idents_grant <- c("Ionocytes", "Esophageal", "Ionocytes", "Goblet cells",
  "Enterocytes", "preEEC", "preSec", "Esophageal", 
  "Enteroendocrine", "Enterocytes",
  "Esophageal", "Enteroendocrine", "Enteroendocrine", "Goblet cells", 
  "Goblet cells", "Ionocytes", "Enterocytes",
  "Enterocytes", "Neutrophils/Macrophages", 
  "ISCs", "Exocrine Pancrease", "Enteroendocrine", "Enteroendocrine",
  "Pharyngeal", "Pharyngeal", "RBCs", "Pharyngeal")

new.idents_grant <- c("IECs", "Esophageal", "IECs", "IECs",
  "IECs", "IECs", "IECs", "Esophageal", 
  "IECs", "IECs",
  "IECs", "IECs", "IECs", "IECs", 
  "IECs", "IECs", "IECs",
  "IECs", "Other", 
  "Other", "Other", "IECs", "IECs",
  "Pharyngeal", "Pharyngeal", "Other", "Pharyngeal")


names(new.idents_grant) <- levels(df_integrated)
df_integrated <- RenameIdents(df_integrated, new.idents_grant)
DimPlot(df_integrated, pt.size = 1, label = TRUE) + scale_color_manual(
  values = colorRampPalette(RColorBrewer::brewer.pal(12, "Dark2"))(length(unique(new.idents_grant)))) +
  theme(legend.position = "none")

ggsave("TAGC-scRNA_2020-04-17.pdf", height = 6, width = 6)

FindAllMarkers(df_integrated) %>% 
  write_csv("scData/WTReferenced_grantMarkers_52PCs35DimIntsres0.8.csv")


# preview a file that would be created by ggsave()
ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

ggpreview(height = 2.5, width = 2.25, units = "in")

library(ggbeeswarm)
ggplot(data = Fulldf3 %>% filter(cluster %in% c("Enterocytes (segment 1)",
  "Enterocytes (segment 2)", "Enterocytes (segment 3)")),
             aes(x = as.factor(cluster), y = MutWt_lfc), size = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 0.5) +
  # stat_summary(geom = "point", fun = median) +
  # geom_beeswarm(data = Fulldf3 %>% filter(p_val_adj < 0.05,
  #                                      avg_logFC > log10(2),
  #                                      gene %in% anterior),
  #            aes(x = as.factor(cluster), y = MutWt_lfc), groupOnX = TRUE,
  #   color = "magenta", size = 0.5) +
  # geom_beeswarm(data = Fulldf3 %>% filter(p_val_adj < 0.05,
  #                                      avg_logFC > log10(2),
  #                                      gene %in% ileal),
  #            aes(x = as.factor(cluster), y = MutWt_lfc), groupOnX = TRUE,
  #   color = "green4", size = 0.5) +
  # geom_beeswarm(data = Fulldf3 %>% filter(p_val_adj < 0.05,
  #                                      avg_logFC > log10(2),
  #                                      gene %in% colon),
  #            aes(x = as.factor(cluster), y = MutWt_lfc), groupOnX = TRUE,
  #   color = "orange", size = 0.5) +
  labs(
    x = "Enterocyte Clusters",
    y = expression(
      log[2] * bgroup("(", frac('uhrf1 Mutant', Wildtype), ")") * " RPKMs")
  ) + scale_x_discrete(label = c("Anterior - Jej", "Anterior - Ileal", "Posterior")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, hjust = 1, angle = 30, color = "black")) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey30")
ggpreview(height = 2.5, width = 2.5, units = "in")
ggsave(
  "C:/Users/drako/Dropbox/Rscripts/F31 FIGURES/Bagnat_scRNASeq_Refined.pdf", height = 2.5, width = 2.5, units = "in")

### End Grant Figures
Fulldf3 %>% filter(p_val_adj < 0.05,
  avg_logFC > log10(2)) %>% group_by(cluster) %>% count() %>% summary()

DimPlot(df_integrated, pt.size = 0.2) + 
  theme(axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "none")


# Saving a matrix of counts per gene per cell
write.table(as.data.frame(as.matrix(df_integrated@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)

df <- read_tsv('Gene_Count_per_Cell_scRNA.tsv')

clust17_0 <- subset(df_integrated.cluster17, idents = '0')
clust17_1 <- subset(df_integrated.cluster17, idents = '1')

clust17_0WT <- subset(clust17_0, subset = orig.ident == 'WT')
clust17_1WT <- subset(clust17_1, subset = orig.ident == 'WT')
clust17_0MUT <- subset(clust17_0, subset = orig.ident == 'MUT')
clust17_1MUT <- subset(clust17_1, subset = orig.ident == 'MUT')



write.table(as.data.frame(as.matrix(clust17_0@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA_17_0.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)
write.table(as.data.frame(as.matrix(clust17_1@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA_17_1.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)

write.table(as.data.frame(as.matrix(clust17_0WT@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA_17_0WT.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)
write.table(as.data.frame(as.matrix(clust17_1WT@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA_17_1WT.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)
write.table(as.data.frame(as.matrix(clust17_0MUT@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA_17_0MUT.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)
write.table(as.data.frame(as.matrix(clust17_1MUT@assays[["RNA"]]@data)) %>%
    rownames_to_column("Gene"),
  file = 'Gene_Count_per_Cell_scRNA_17_1MUT.tsv',
  quote = FALSE, sep = '\t', col.names = NA, row.names = TRUE)


refMarkers <- read_csv("scData/WTReferenced_Markers_52PCs35DimIntsres0.8.csv")
predGenes <- c("syt7b", "fer1l6", "itga2.1", "fabp1b", "fabp2", "ada", "rbp2a", "chia.2",
  "apoa1b", "apoa1a", "pyyb", "ccka", "penka", "adcyap1a", "gcga", "galn",
  "xnpep2", "ghrl", "sst1.2", "sst2", "pnzr1", "si:ch73-14h1.2", "malb", "si:zfos-2372e4.1",
  "s100a1", "FP236331.1", "pax4", "urahb", "bdnf", "trpa1b", "si:dkey-234i14.2",
  "lmx1ba", "hbegfb", "tph1b", "si:dkey-61f9.1", "tnfrsf11b", "adgrfb3a", "pou2f3",
  "si:dkey-203a12.9", "edrf1", "hsp70.1", "fabp6", "slc10a2", "spi1", "twist1a",
  "snai1a", "pmp22a", "ela3l", "pdia2", "prss59.2", "prss1", "cbp1", "amy2a",
  "CR556712.1", "calca", "ddhd2", "insl5a", "gcga", "insl5b", "tac3a", "galn", "btc",
  "vipb", "cldni", "cxl34b.11", "cx30.3", "si:dkey-96g2.1", "shroom4", "edn2", "tnfb",
  "dhrs13a.2", "cldne", "anxa1c")


## Looking at Transcription factor binding sites and average gene expression per cluster-----
# This is from the TF database
DanioTFs <- read_tsv("Danio_rerio_TF.txt") %>% rename("gene" = "Symbol")
dfExpression <- read_csv("AverageGeneExpression_FinalCluster.csv") %>% rename("gene" = "X1")
clust17_1UP <- read_tsv("knownResults-homer_Clust17-1UP.txt")
clust17_1UP$`Motif Name` <- clust17_1UP$`Motif Name` %>%
  str_remove("\\(.*") %>% str_to_lower() # Removes the excess stuff in Homer output (may want that though?)

sum(clust17_1UP$`Motif Name` %in% dfExpression$gene)

ggplot(filter(dfExpression, gene %in% DanioTFs$gene)) +
  geom_point(aes(x = RNA.17, y = RNA.0))

dfExpressionLong <- dfExpression %>% pivot_longer(cols = 2:28)
dfExpressionLongSummary <- dfExpressionLong %>% 
  group_by(gene) %>% 
  summarize(median = median(value), sd = stats::sd(value))
dfExpressionLong <- inner_join(dfExpressionLong, dfExpressionLongSummary) %>% 
  mutate(zscore = (value - median)/sd)

dfExpressionLong <- dfExpressionLong %>% group_by(gene) %>% filter(median != 0)

ggplot(filter(dfExpressionLong, gene %in% DanioTFs$gene)) +
  geom_tile(aes(x = name, fill = zscore, y = gene))

dfMatrix <- dfExpressionLong %>% select(gene, name, value)
dfMatrix <- dfMatrix %>% pivot_wider(names_from = name, values_from = value)
dfMatrix <- column_to_rownames(dfMatrix, var = "gene") %>%
  as.matrix()
dfMatrix <- scale(t(dfMatrix))

d <- dist(dfMatrix)
hc <- hclust(d)
plot(hc)
library(pheatmap)
pheatmap(dfMatrix, show_colnames = FALSE)

heatmap(dfMatrix, scale = "none")



# Using the Adult zebrafish FAIRE-Seq with zCNEs 10kb upstream and downstream of genes
library(biomaRt)

ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
datasets
searchDatasets(mart = ensembl, pattern = "drerio")
ensembl <- useDataset("drerio_gene_ensembl", mart = ensembl)
listAttributes(ensembl)[1:50,]
listFilters(ensembl)[1:10,]

df <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id',
  'ensemble_transcript_id_version', 'entrezgene_id'), 
        filters = 'name', 
        mart = ensembl)
ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
#for danRer10
ensembl91 <- useMart(host = 'dec2017.archive.ensembl.org', 
                     biomart = 'ensembl', 
                     dataset = 'drerio_gene_ensembl')
head(listAttributes(ensembl91), 10)
searchAttributes(ensembl, pattern = "GO")

bm_list <- getBM(
  c("ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "ensembl_transcript_id_version",
    "zfin_id_symbol"),
  mart = ensembl91)

bm_list <- getBM(
  c("ensembl_gene_id",
    "external_gene_name",
    "go_id",
    "name_1006"),
  mart = ensembl)

peaks <- read_tsv("SelectedPeaks_Final_toMerge.txt",
  col_names = c("Chromosome", "Start", "End", "Gene", "NA", "Strand"))
peaksUp <- read_tsv("UPstream_Gene.bed",
  col_names = c("Chromosome", "Start", "End", "Gene", "NA", "Strand"))
peaksUp <- peaksUp %>% distinct()


peaksDown <- read_tsv("DownStream.bed",
  col_names = c("Chromosome", "Start", "End", "Gene", "NA", "Strand"))

peaksDown$Gene <- str_replace(peaksDown$Gene, "_[:graph:]+", "")
peaksDown <- peaksDown %>% distinct()

AllPeaks <- read_tsv("10kbUPDOWN_Gene.bed",
  col_names = c("Chromosome", "Start", "End", "Gene", "NA", "Strand"))

AllPeaks$Gene <- str_replace(AllPeaks$Gene, "_[:graph:]+", "")
AllPeaks <- AllPeaks %>% distinct()


# This is the final file
# Combining FAIRE-Seq data with TF binding motifs annotations----
FinalFairePeaks <- read_tsv("FAIRE-Seq_IDA/FinalFAIRE_10kbwindows_4selection.bed",
  col_names = c("Chromosome", "Start", "End", "Peak",
    "ChromosomeG", "StartG", "EndG", "ensembl_transcript_id_version",
    "NA", "Strand"))
FinalFairePeaksAnno <- FinalFairePeaks %>%
  left_join(bm_list91 %>% dplyr::select(ensembl_gene_id,
    external_gene_name, 
    ensembl_transcript_id_version)) %>% dplyr::select(1:4, 11, 12) %>% distinct()

# Query the average interval for these FAIRE-Peaks
library(GenomicRanges)
simpleRanges <- GRanges(seqnames = FinalFairePeaksAnno$Chromosome,
  ranges = IRanges(start = FinalFairePeaksAnno$Start,
    end = FinalFairePeaksAnno$End),
  names = FinalFairePeaksAnno$Peak)

FinalFairePeaksAnno <- FinalFairePeaksAnno %>%
  mutate(size = End - Start) %>% filter(size > 20)

FinalFairePeaksAnno130 <- FinalFairePeaksAnno %>%
  mutate(size = End - Start) %>% filter(size > 130)


Up16 <- read_csv("Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_16_FinalUMAP.csv") %>% 
  filter(avg_logFC > 0)
Dw16 <- read_csv("Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_16_FinalUMAP.csv") %>% 
  filter(avg_logFC < 0)

UpPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Up16$gene)
DownPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Dw16$gene)

write_tsv(UpPeaks %>% dplyr::select(1:3,6), "16_10kb_CNE-FAIRE_Upreg_fixed_names.bed", col_names = FALSE)
write_tsv(DownPeaks %>% dplyr::select(1:3,6), "16_10kb_CNE-FAIRE_Downreg_fixed_names.bed", col_names = FALSE)


Up4 <- read_csv("Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_4_FinalUMAP.csv") %>% 
  filter(avg_logFC > 0)
Dw4 <- read_csv("Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_4_FinalUMAP.csv") %>% 
  filter(avg_logFC < 0)

UpPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Up4$gene)
DownPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Dw4$gene)

write_tsv(UpPeaks %>% dplyr::select(1:3), "4_10kb_CNE-FAIRE_Upreg_fixed_names.bed", col_names = FALSE)
write_tsv(DownPeaks %>% dplyr::select(1:3), "4_10kb_CNE-FAIRE_Downreg_fixed_names.bed", col_names = FALSE)


DownPEnterocytes <- c(Dw4$gene, Dw16$gene)
UpPEnterocytes <- c(Up4$gene, Up16$gene)

UpPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% UpPEnterocytes)
DownPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% DownPEnterocytes)

write_tsv(UpPeaks %>% dplyr::select(1:3), "Enterocytes_10kb_CNE-FAIRE_Upreg_fixed.bed", col_names = FALSE)
write_tsv(DownPeaks %>% dplyr::select(1:3), "Enterocytes_10kb_CNE-FAIRE_Downreg_fixed.bed", col_names = FALSE)


Up17 <- read_csv("Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_17_FinalUMAP.csv") %>% 
  filter(avg_logFC > 0)
Dw17 <- read_csv("Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_17_FinalUMAP.csv") %>% 
  filter(avg_logFC < 0)

UpPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Up17$gene)
DownPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Dw17$gene)

write_tsv(UpPeaks %>% dplyr::select(1:3, 6), "17_10kb_CNE-FAIRE_Upreg_fixed_names.bed", col_names = FALSE)
write_tsv(DownPeaks %>% dplyr::select(1:3,6), "17_10kb_CNE-FAIRE_Downreg_fixed_names.bed", col_names = FALSE)

# write_tsv(peaksUp, "10kb_Up_w-Gene.bed", col_names = FALSE)
# write_tsv(peaksDown, "10kb_Down_no-Gene.bed", col_names = FALSE)
# 
# write_tsv(AllPeaks, "10kb_UpDown_w-Gene_cat.bed", col_names = FALSE)

temporaryUp <- read_tsv("Enterocytes_10kb_CNE-FAIRE_Upreg_fixed.bed")
temporaryDown <- read_tsv("Enterocytes_10kb_CNE-FAIRE_Downreg_fixed.bed")

Up17_1 <- read_csv("Clust17_2subclusterRes0.5_1_MUTvWTdiff.csv") %>% 
  filter(avg_logFC > 0)
Dw17_1 <- read_csv("Clust17_2subclusterRes0.5_1_MUTvWTdiff.csv") %>% 
  filter(avg_logFC < 0)


UpPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Up17_1$gene)
DownPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Dw17_1$gene)

write_tsv(UpPeaks %>% dplyr::select(1:3), "17_1_10kb_CNE-FAIRE_Upreg_fixed.bed", col_names = FALSE)
write_tsv(DownPeaks %>% dplyr::select(1:3), "17_1_10kb_CNE-FAIRE_Downreg_fixed.bed", col_names = FALSE)




Up17_0 <- read_csv("Clust17_2subclusterRes0.5_0_MUTvWTdiff.csv") %>% 
  filter(avg_logFC > 0)
Dw17_0 <- read_csv("Clust17_2subclusterRes0.5_0_MUTvWTdiff.csv") %>% 
  filter(avg_logFC < 0)


UpPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Up17_0$gene)
DownPeaks <- FinalFairePeaksAnno %>% filter(external_gene_name %in% Dw17_0$gene)

write_tsv(UpPeaks %>% dplyr::select(1:3), "17_0_10kb_CNE-FAIRE_Upreg_fixed.bed", col_names = FALSE)
write_tsv(DownPeaks %>% dplyr::select(1:3), "17_0_10kb_CNE-FAIRE_Downreg_fixed.bed", col_names = FALSE)



makeBEDcluster <- function(cluster) {

  cluster <- as.character(cluster)
  
  UpCluster <- read_csv("MutvWT/MutvWT_Cluster_", cluster, "_FinalUMAP.csv") %>% 
    filter(avg_logFC > 0)
  
  DwCluster <- read_csv("MutvWT/MutvWT_Cluster_", cluster, "_FinalUMAP.csv") %>% 
    filter(avg_logFC < 0)

  UpPeaks <- FinalFairePeaksAnno %>% 
    filter(external_gene_name %in%
        UpCluster$gene)
  
  DownPeaks <- FinalFairePeaksAnno %>% 
    filter(external_gene_name %in%
        DwCluster$gene)

  write_tsv(UpPeaks %>% dplyr::select(1:3),
    paste0(cluster, "_10kb_CNE-FAIRE_Upreg_fixed.bed"), col_names = FALSE)
  write_tsv(DownPeaks %>% dplyr::select(1:3),
    paste0(cluster, "_10kb_CNE-FAIRE_Downreg_fixed.bed"), col_names = FALSE)

}

# Prepare the random background file
FinalBkgFairePeaks <- read_tsv("FinalFAIRE_10kbwindows_4background.bed",
  col_names = c("Chromosome", "Start", "End",
    "ChromosomeG", "StartG", "EndG", "ensembl_transcript_id_version"))
FinalBkgFairePeaksAnno <- FinalBkgFairePeaks %>%
  left_join(bm_list %>% dplyr::select(ensembl_gene_id,
    external_gene_name, 
    ensembl_transcript_id_version)) %>% dplyr::select(1:3, 8,9) %>% distinct()

set.seed(54321)
BckgPeaks <- FinalBkgFairePeaksAnno[sample(nrow(FinalBkgFairePeaksAnno), 800), ]

write_tsv(BckgPeaks %>% dplyr::select(1:3), "17_10kb_CNE-FAIRE_DownRegBackground.bed", col_names = FALSE)
write_tsv(BckgPeaks %>% dplyr::select(1:3), "17_10kb_CNE-FAIRE_UpRegBackground.bed", col_names = FALSE)



DefaultAssay(df_integrated.cluster17) <- "RNA"
# Single cell heatmap of feature expression
DoHeatmap(df_integrated, features = predGenes, size = 3)

FeaturePlot(df_integrated, features = c("fabp6", "slc10a2"), min.cutoff = "q10", max.cutoff = "q90", split.by = "library")

clust17_avg <- AverageExpression(df_integrated.cluster17, assays = "RNA")
write.csv(clust17_avg, "Cluster17_subclusterAverageExpression.csv")


cNumbersAll <- as.data.frame(table(Idents(df_integrated))) %>%
  mutate(FreqWT = as.data.frame(table(Idents(subset(df_integrated, 
    subset = orig.ident == "WT"))))[[2]],
    FreqMUT = as.data.frame(table(Idents(subset(df_integrated, 
    subset = orig.ident == "MUT"))))[[2]])

library(gt)

dgt <- gt(data = cNumbersAll) %>%
  tab_header(
    title = "Cell Numbers by Cluster") %>%
  cols_label(Var1 = "Cluster", Freq = "Total", FreqWT = "Wild-type", FreqMUT = "Mutant")
gtsave(dgt,"CellNumberTable.pdf")


clust17_0diff <- read_csv("Clust17_2subclusterRes0.5_0_MUTvWTdiff.csv")
clust17_1diff <- read_csv("Clust17_2subclusterRes0.5_1_MUTvWTdiff.csv")
tibble(summarise(clust17_0diff,
  Down_0 = sum(avg_logFC < 0), Up_0 = sum(avg_logFC > 0)),
summarise(clust17_1diff,
  Down_1 = sum(avg_logFC < 0), Up_1 = sum(avg_logFC > 0))
)

tempdf <- tibble(
 Subcluster = c(0, 1, 0, 1),
  Direction = c("Up", "Up", "Down", "Down"),
  Count = c(355, 1030, 299, 477)
)


ggplot(tempdf, aes(x = Subcluster, y = Count)) +
  geom_col(aes(fill = Direction)) + 
  scale_fill_manual(values = c("grey30", "grey70")) +
  theme_bw() + ggtitle("Absolute Differential Gene Count in Cluster 17")


tempdf <- tempdf %>% group_by(Subcluster) %>%
  mutate(RelCount = Count/sum(Count))

ggplot(tempdf, aes(x = Subcluster, y = RelCount)) +
  geom_col(aes(fill = Direction)) + 
  scale_fill_manual(values = c("grey30", "grey70")) +
  theme_bw() + ggtitle("Relative Differential Gene Count in Cluster 17")


# Clustering Clusters by average expression and differential mut/wt expression
avg_exp <- read_csv("AverageGeneExpression_FinalCluster.csv")
df_mvw <- read_csv("MutvWT/MutvWT_AllClusters_FinalUMAP.csv")

# Maybe select variable features first (lots of zeros here)
features <- row.names(df_integrated@reductions$pca@feature.loadings)
features <- bigUniqueFindTest2$gene

# Make df_mvw a wide matrix
df_mvw <- df_mvw %>% select(gene, avg_logFC, Cluster) %>% 
  pivot_wider(id_cols = gene, names_from = Cluster, values_from = avg_logFC)
df_mvw <- df_mvw %>% replace(is.na(.), 0)
df_mvw <- df_mvw %>% column_to_rownames("gene") %>% as.matrix()

# Also make avg_exp a matrix
avg_exp <- avg_exp %>% filter(X1 %in% features) %>% column_to_rownames("X1") %>% as.matrix()

dmat_express <- scale(avg_exp)
heatmap(dmat_express, scale = "none")
heatmap(avg_exp, scale = "none")
heatmap(avg_exp, scale = "row") # Row-wise normalization (where are certain genes expressed highest?)
heatmap(avg_exp, scale = "column") # Column-wise normalization (profile of expression within each cluster)

heatmap(df_mvw, scale = "none")
heatmap(df_mvw, scale = "row")

cormat <- round(cor(avg_exp),2)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
ggplot(data = melted_cormat, aes(x = as.factor(Var1), y = as.factor(Var2), fill = value)) + 
  geom_tile() + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Function!!!
reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1 - cormat)/2)
hc <- hclust(dd)
cormat <- cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)

plot(hclust(dist(t(df_mvw))))
df_mvwClustered <- hclust(dist(df_mvw))

plot(hclust(dist(t(avg_exp))))
plot(hclust(dist(avg_exp)))

# Function to write all FindMarkerComparisons into a folder----
FindAllUniqueMarkers <- function(object, nclusters) {
  for (i in 0:nclusters) {
    j <- i + 1
    while (j <= nclusters) {
        print(paste0("Comparing clusters ", i, " & ", j))
        temp <- FindMarkers(object, ident.1 = i, ident.2 = j, min.pct.diff = 0.25) %>%
          rownames_to_column(var = "gene")
        temp %>% filter(avg_logFC < 0) %>%
          mutate(avg_logFC = -1 * avg_logFC,
            cluster1 = j, cluster2 = i) %>%
          write_csv(paste0("UniqueMarkers2/FindMarkers_", j, "_vs_", i, "_minpctdiff0.25.csv"))
        temp %>% filter(avg_logFC > 0) %>%
          mutate(avg_logFC = avg_logFC,
            cluster1 = i, cluster2 = j) %>%
          write_csv(paste0("UniqueMarkers2/FindMarkers_", i, "_vs_", j, "_minpctdiff0.25.csv"))
        j <- j + 1
    }
  }
}

# Function to compare all clusters to superclusters----
FindAllSuperMarkers <- function(object, nclusters) {
  for (i in 0:nclusters) {
    # Ionocytes
        print(paste0("Comparing clusters ", i, " & Ionocytes"))
        temp <- FindMarkers(object, ident.1 = i, ident.2 = c(0, 2, 15)) %>%
          rownames_to_column(var = "gene")
        temp %>% filter(avg_logFC < 0) %>%
          mutate(avg_logFC = -1 * avg_logFC,
            cluster1 = "Ionocytes", cluster2 = i) %>%
          write_csv(paste0("SuperClusters/FindMarkers_Ionocytes_vs_", i, ".csv"))
        temp %>% filter(avg_logFC > 0) %>%
          mutate(avg_logFC = avg_logFC,
            cluster1 = i, cluster2 = "Ionocytes") %>%
          write_csv(paste0("SuperClusters/FindMarkers_", i, "_vs_Ionocytes.csv"))
    # Enterocytes
        print(paste0("Comparing clusters ", i, " & Enterocytes"))
        temp <- FindMarkers(object, ident.1 = i, ident.2 = c(4, 16, 17, 9)) %>%
          rownames_to_column(var = "gene")
        temp %>% filter(avg_logFC < 0) %>%
          mutate(avg_logFC = -1 * avg_logFC,
            cluster1 = "Enterocytes", cluster2 = i) %>%
          write_csv(paste0("SuperClusters/FindMarkers_Enterocytes_vs_", i, ".csv"))
        temp %>% filter(avg_logFC > 0) %>%
          mutate(avg_logFC = avg_logFC,
            cluster1 = i, cluster2 = "Enterocytes") %>%
          write_csv(paste0("SuperClusters/FindMarkers_", i, "_vs_Enterocytes.csv"))
    # Goblet Cells
        print(paste0("Comparing clusters ", i, " & GobletCells"))
        temp <- FindMarkers(object, ident.1 = i, ident.2 = c(3, 13, 14)) %>%
          rownames_to_column(var = "gene")
        temp %>% filter(avg_logFC < 0) %>%
          mutate(avg_logFC = -1 * avg_logFC,
            cluster1 = "GobletCells", cluster2 = i) %>%
          write_csv(paste0("SuperClusters/FindMarkers_GobletCells_vs_", i, ".csv"))
        temp %>% filter(avg_logFC > 0) %>%
          mutate(avg_logFC = avg_logFC,
            cluster1 = i, cluster2 = "GobletCells") %>%
          write_csv(paste0("SuperClusters/FindMarkers_", i, "_vs_GobletCells.csv"))
    # EECs
        print(paste0("Comparing clusters ", i, " & EECs"))
        temp <- FindMarkers(object, ident.1 = i, ident.2 = c(8, 11, 12, 21, 22)) %>%
          rownames_to_column(var = "gene")
        temp %>% filter(avg_logFC < 0) %>%
          mutate(avg_logFC = -1 * avg_logFC,
            cluster1 = "EECs", cluster2 = i) %>%
          write_csv(paste0("SuperClusters/FindMarkers_EECs_vs_", i, ".csv"))
        temp %>% filter(avg_logFC > 0) %>%
          mutate(avg_logFC = avg_logFC,
            cluster1 = i, cluster2 = "EECs") %>%
          write_csv(paste0("SuperClusters/FindMarkers_", i, "_vs_EECs.csv"))
    # Esophageal
        print(paste0("Comparing clusters ", i, " & Esophageal"))
        temp <- FindMarkers(object, ident.1 = i, ident.2 = c(1, 7, 10)) %>%
          rownames_to_column(var = "gene")
        temp %>% filter(avg_logFC < 0) %>%
          mutate(avg_logFC = -1 * avg_logFC,
            cluster1 = "Esophageal", cluster2 = i) %>%
          write_csv(paste0("SuperClusters/FindMarkers_Esophageal_vs_", i, ".csv"))
        temp %>% filter(avg_logFC > 0) %>%
          mutate(avg_logFC = avg_logFC,
            cluster1 = i, cluster2 = "Esophageal") %>%
          write_csv(paste0("SuperClusters/FindMarkers_", i, "_vs_Esophageal.csv"))
  }
}



# Test the loops
  for (i in 0:3) {
    for (j in 0:3) {
      if (i != j) {
        print(paste0("Comparing clusters ", i, " & ", j))
      } else {
        print("Redundant cluster comparison, skipping")
      }
    }
  }

# Do the run on the real data
FindAllUniqueMarkers(df_integrated, nclusters = 26)

# Read in all files
library(fs)
csv_files <- fs::dir_ls("UniqueMarkers")
csv_files <- fs::dir_ls("UniqueMarkers2")
csv_files <- fs::dir_ls("SuperClusters")


read_Cnames <- function(filename) {
  read_csv(filename, col_types = cols(cluster1 = col_character(), cluster2 = col_character()))
}

bigUniqueFind <- csv_files %>% 
  map_dfr(read_Cnames, .id = "source")

bigUniqueFind <- csv_files %>% 
  map_dfr(read_csv, .id = "source")

bigUniqueFindTest <- bigUniqueFind %>% 
   select(gene, cluster1) %>% distinct() %>%
  count(gene) %>% filter(n == 1) %>% left_join(bigUniqueFind %>%
  filter(avg_logFC > 0) %>% select(gene, cluster1)) %>% distinct() %>%
  mutate(cluster1 = as.numeric(cluster1)) %>%
  arrange(cluster1)

bigUniqueFindTest2 <- bigUniqueFind %>% mutate(diff.pct = pct.1 - pct.2) %>%
  filter(diff.pct > 0.1) %>%
  select(gene, cluster1) %>% count(gene, cluster1)  %>%
  filter(n == 26) %>%
  mutate(cluster1 = as.numeric(cluster1)) %>%
  arrange(cluster1)
  
bigUniqueFindTest2 <- bigUniqueFind %>% 
  select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(4, 16, 9, 17),
    !(cluster2 %in% c(4, 16, 9, 17))) %>% 
  count(gene, cluster1) %>% filter(n == 23) %>%
  count(gene, name = "totalTimes")

bigUniqueFindTest2 <- bigUniqueFind %>% 
  select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(4, 16),
    !(cluster2 %in% c(4, 16))) %>% 
  count(gene, cluster1) %>% filter(n == 25) %>%
  count(gene, name = "totalTimes") %>% filter(totalTimes == 2)


bigUniqueFindTest2 <- bigUniqueFind %>% 
  select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(9, 17),
    !(cluster2 %in% c(9, 17))) %>% 
  count(gene, cluster1) %>% filter(n == 23) %>%
  count(gene, name = "totalTimes") %>% filter(totalTimes == 2)

bigUniqueFindTest2 <- bigUniqueFind %>% 
  select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(4, 16, 9),
    !(cluster2 %in% c(4, 16, 9))) %>% 
  count(gene, cluster1) %>% filter(n == 24) %>%
  count(gene, name = "totalTimes") %>% filter(totalTimes == 3)

bigUniqueFindTest2 <- bigUniqueFind %>% 
  select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(4, 16, 17),
    !(cluster2 %in% c(4, 16, 17))) %>% 
  count(gene, cluster1) %>% filter(n == 24) %>%
  count(gene, name = "totalTimes") %>% filter(totalTimes == 3)


hist(as.numeric(bigUniqueFindTest2$n))
table(bigUniqueFindTest2$cluster1)

write_csv(bigUniqueFindTest2, "UniqueMarkers/AllUniqueMarkers.csv")
bigUniqueFindTest2 <- read_csv("UniqueMarkers/AllUniqueMarkers.csv")

bigUniqueFind %>% select(gene, cluster1) %>% filter(gene == "abat")

for (i in 1:nrow(bigUniqueFindTest2)) {
  print(paste("Plotting Unique Gene: ",
    bigUniqueFindTest2$gene[i],
    " for Cluster ",
    bigUniqueFindTest2$cluster1[i],
    ", ", i, "/", nrow(bigUniqueFindTest2)))
FeaturePlot(df_integrated,
  features = bigUniqueFindTest2$gene[i],
  label = TRUE) + 
  ggtitle(label = paste0(bigUniqueFindTest2$gene[i],
    ", Marker for Cluster ",
    bigUniqueFindTest2$cluster1[i]) )
gene_name <- as.character(bigUniqueFindTest2$gene[i]) %>%
    str_replace(":", "_")
ggsave(paste0("UniqueMarkersFeaturePlots/FeaturePlot_",
  gene_name,
  "_forCluster_",
  bigUniqueFindTest2$cluster1[i],
  ".pdf"), height = 6, width = 6)
}


FindAllMarkers(df_integrated, min.pct =  0.25, 
                             min.diff.pct = 0.25, only.pos = TRUE) %>%
  write.csv("FindAllMarkers.Wilcox_mindiffpct_0.5.csv")

nrow(mtcars)
mtcars %>% glimpse()


df_sample <- data.frame(Cluster = 0:26, Value = c(15, rep(20:25, 4), 15, 20))

smalldftest <- list()
for (i in 0:(nrow(df_sample) - 1)) {
  j <- i + 1
  while (j <= nrow(df_sample)) {
    print(paste0(i,"_",j))
    smalldftest[[nrow(df_sample)*i + j]] <- paste0(i,"_",j)
    j <- j + 1
  }
}
smalldftest 
df_sample %>% mutate(Value = -1*Value)

# Addressing some scientific questions during edits-----
# Question 3
DefaultAssay(df_integrated.cluster17) <- "RNA"

FeaturePlot(df_integrated.cluster17,
  features = c("fabp6", "amn"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_fabp6-amn.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("fabp6", "cubn"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_fabp6-cubn.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("slc10a2", "amn"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_slc10a2-amn.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("slc10a2", "cubn"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_slc10a2-cubn.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("fabp6", "dab2"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_fabp6-dab2.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("slc10a2", "dab2"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_slc10a2-dab2.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("fabp6", "pllp"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_fabp6-pllp.pdf", height = 5, width = 20)

FeaturePlot(df_integrated.cluster17,
  features = c("slc10a2", "pllp"), blend = TRUE, cols = c("grey80", "red", "blue"),
  blend.threshold = 0.1)
ggsave("BlendedFeaturePlot_slc10a2-pllp.pdf", height = 5, width = 20)



# Feature Scatter
FeatureScatter(df_integrated.cluster17,
  feature1 = "fabp6", feature2 = "amn")
ggsave("FeatureScatter_fabp6-amn.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "fabp6", feature2 = "cubn")
ggsave("BlendedFeatureScatter_fabp6-cubn.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "slc10a2", feature2 ="amn")
ggsave("FeatureScatter_slc10a2-amn.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "slc10a2", feature2 = "cubn")
ggsave("FeatureScatter_slc10a2-cubn.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "fabp6", feature2 = "dab2")
ggsave("BlendedFeatureScatter_fabp6-dab2.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "slc10a2", feature2 = "dab2")
ggsave("FeatureScatter_slc10a2-dab2.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "fabp6", feature2 = "pllp")
ggsave("BlendedFeatureScatter_fabp6-pllp.pdf", height = 5, width = 20)

FeatureScatter(df_integrated.cluster17,
  feature1 = "slc10a2", feature2 = "pllp")
ggsave("FeatureScatter_slc10a2-pllp.pdf", height = 5, width = 20)


##
## Subsetting cluster 12
##
df_integrated.cluster12 <- subset(df_integrated, idents = '12')
DefaultAssay(df_integrated.cluster12) <- "integrated"
df_integrated.cluster12 <- FindVariableFeatures(df_integrated.cluster12,
                                                verbose = FALSE)
df_integrated.cluster12 <- ScaleData(df_integrated.cluster12, verbose = FALSE)
df_integrated.cluster12 <- RunPCA(object = df_integrated.cluster12,  npcs = 100)
df_integrated.cluster12 <- ProjectDim(object = df_integrated.cluster12)

df_integrated.cluster12 <- JackStraw(object = df_integrated.cluster12,
                                     num.replicate = 100,
                                     dims = 100,
                                     prop.freq = 0.1, verbose = TRUE)
df_integrated.cluster12 <- ScoreJackStraw(object = df_integrated.cluster12,
                                          dims = 1:100) 
JackStrawPlot(object = df_integrated.cluster12, dims = 1:12)

pcs.use = 8

df_integrated.cluster12 <- FindNeighbors(object = df_integrated.cluster12,
                                         dims = 1:pcs.use)
df_integrated.cluster12 <- FindClusters(object = df_integrated.cluster12,
                                        resolution = 0.5)
df_integrated.cluster12 <- RunUMAP(object = df_integrated.cluster12,
                                   dims = 1:pcs.use,
                                   assay = "integrated")
DimPlot(df_integrated.cluster12) + ggtitle(label = "Sub-clustering Cluster 12")


# Conserved Markers
get_conserved <- function(cluster){
  FindConservedMarkers(df_integrated,
                       ident.1 = cluster,
                       grouping.var = "library",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(bm_list[, c("external_gene_name", "description")]),
               by = c("gene" = "external_gene_name")) %>%
    cbind(cluster_id = cluster, .)
  }

conserved_markers <- map_dfr(0:26, get_conserved)
conserved_markers <- conserved_markers %>%
  mutate(pct_diff = ((WT_pct.1 + MUT_pct.1) - (WT_pct.2 + MUT_pct.2))/2,
    logFC_overall = 0.5*(WT_avg_logFC + MUT_avg_logFC))

conserved_markers10 <- conserved_markers %>% group_by(cluster_id) %>%
  top_n(10, wt = logFC_overall) %>% arrange(cluster_id, desc(pct_diff), logFC_overall)


# Find Marker Genes
ccyle <- c("ccna2", "ccnb1", "ccnb1", "ccne",
  "p21", "rb")

FindAllMarkers(df_integrated,
  test.use = "negbinom",
  min.pct = 0.2,
  only.pos = TRUE) %>%
  write.csv("FindAllMarkers_negBinomial_min.pct_0.2.csv")




print(df_integrated[["pca"]], dims = 10, nfeatures = 10)

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:25),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(df_integrated, 
                     vars = columns)

umap_label <- FetchData(df_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x = mean(UMAP_1), y = mean(UMAP_2))
  
# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:25), function(pc){
        ggplot(pc_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color = pc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data = umap_label, 
                          aes(label = ident, x, y)) +
                ggtitle(pc)
}) %>% 
        cowplot::plot_grid(plotlist = .)



# Do Haber dataset
library(Matrix)

zz = gzfile("GSE92332_atlas_UMIcounts.txt.gz",'rt')  
raw_counts <- read.table(zz, 
  sep = "\t")

head(raw_counts)
mIEC <- CreateSeuratObject(raw_counts,
  project = "mouseIEC_scRNAseq")


mIEC$cell.type <- str_extract(mIEC@assays$RNA@counts@Dimnames[[2]],
  "([:alnum:]|\\.)+$")

mIEC <- subset(mIEC, subset = nCount_RNA < 30000 & nFeature_RNA < 3000)

mIEC <- NormalizeData(mIEC, verbose = FALSE)
mIEC <- FindVariableFeatures(mIEC, verbose = FALSE)
mIEC <- ScaleData(object = mIEC, verbose = FALSE)
mIEC <- RunPCA(object = mIEC,
  features = VariableFeatures(object = mIEC),
  npcs = 100)


ElbowPlot(mIEC, ndims = 30)

# Redo Jackstraw for the Haber data
mIEC <- JackStraw(object = mIEC,
                           num.replicate = 100,
                           dims = 100,
                           prop.freq = 0.1,
                           verbose = TRUE)
mIEC <- ScoreJackStraw(object = mIEC,
                                dims = 1:100) 
#
JackStrawPlot(object = mIEC, dims = 1:40)


pcs.use <- 28 # The selected number of PCs

mIEC <- RunUMAP(mIEC,
  reduction = "pca", dims = 1:pcs.use)

mIEC <- FindNeighbors(mIEC,
  reduction = "pca", dims = 1:pcs.use)

mIEC <- FindClusters(mIEC, resolution = 0.65)

DimPlot(object = mIEC, reduction = "umap",
        label = TRUE)

DimPlot(object = mIEC, reduction = "umap",  group.by = "cell.type",
        label = TRUE, label.size = 2)

FindAllMarkers(mIEC, min.pct = 0.2) %>% write.csv("Haber_15Cluster_Markers.csv")

Idents(object = mIEC) <- 'cell.type'

FindAllMarkers(mIEC, min.pct = 0.2) %>% write.csv("Haber_15CellType_Markers.csv")

haber_celltypes <- read_csv("Haber_15CellType_Markers.csv") %>%
  filter(avg_logFC > 0)

haber_celltypes20 <- haber_celltypes %>% group_by(cluster) %>% arrange(cluster, desc(avg_logFC)) %>%
  top_n(20)


AverageExpression(mIEC) %>% write.csv("AverageExpression_Haber_15CellType.csv")

saveRDS(mIEC, "Haber.rds")
mIEC <- readRDS("Haber.rds")

mIEC_exp <- read_csv("AverageExpression_Haber_15CellType.csv")

mIEC_exp <- mIEC_exp %>% column_to_rownames("X1") %>% as.matrix()

heatmap(mIEC_exp, scale = "row")

keeping <- str_match(levels(Idents(mIEC)), "(Enterocyte|TA)[:graph:]+")[,1] %>% as.character() %>%
  na.omit()

mIEC_enterocytes <- subset(mIEC, idents = keeping)
plot <- DimPlot(mIEC_enterocytes)
mIEC_enterocytes <- CellSelector(plot = plot, object = mIEC_enterocytes, ident = 'SelectedCells')

bm_list <-  bm_list %>% filter(mmusculus_homolog_orthology_type == "ortholog_one2one")
haber_celltypes1to1 <- haber_celltypes %>%
  filter(gene %in% bm_list$mmusculus_homolog_associated_gene_name)

haber_celltypes1to1 <- left_join(haber_celltypes1to1, bm_list,
  by = c("gene" = "mmusculus_homolog_associated_gene_name"))

haber_celltypes50 <- haber_celltypes1to1 %>% 
  group_by(cluster) %>% top_n(50, wt = avg_logFC) %>%
  arrange(cluster, desc(avg_logFC)) %>%
  dplyr::select(-X1) %>% write_csv("Haber_1to1ortho_FindMarkers.csv")

# Jia requests 2020-05-28
# (1) 0, 2, 15 shared genes as compared against every other cluster
# (2) differential genes of the super cloud 0/2/15 against everything else.

# (1)
library(fs)
csv_files <- fs::dir_ls("Output/DiffExpTests/Clusters/UniqueMarkers")
csv_files <- fs::dir_ls("UniqueMarkers2")

FAmarkers <- read_csv("Output/DiffExpTests/WTReferenced_Markers_52PCs35DimIntsres0.8.csv")


bigUniqueFind <- csv_files %>% 
  map_dfr(read_csv, .id = "source")

  
bigUniqueFindTest2 <- bigUniqueFind %>% 
  dplyr::select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(0, 2, 15),
    !(cluster2 %in% c(0, 2, 15))) %>% 
  count(gene, cluster1) %>% filter(n == 24) %>%
  count(gene, name = "totalTimes") %>% filter(totalTimes == 3) %>%
  dplyr::select(gene) %>% left_join(bigUniqueFind) %>%
  filter(cluster1 %in% c(0, 2, 15),
    !(cluster2 %in% c(0, 2, 15))) %>% arrange(cluster1, p_val_adj, pct.1) %>%
  dplyr::select(-source)

write_csv(bigUniqueFindTest2, "Unique-0-2-15_FullData.csv")
write_csv(as.data.frame(unique(bigUniqueFindTest2$gene)), "Unique-0-2-15_GeneList.csv")
unique(bigUniqueFindTest2$gene)
# (2)
FindMarkers(df_integrated, ident.1 = c(0, 2, 15), min.pct = 0.2) %>%
  write.csv("FindMarkers_0-2-15_minpct-0.2.csv")


bigUniqueFindTest2 <- bigUniqueFind %>% 
  dplyr::select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(0, 2, 15),
    !(cluster2 %in% c(0, 2, 15))) %>% 
  count(gene, cluster1) %>% filter(n == max(n)) %>%
  dplyr::select(gene) %>% left_join(bigUniqueFind) %>%
  filter(cluster1 %in% c(0, 2, 15)) %>%
  arrange(cluster1, p_val_adj, pct.1) %>%
  dplyr::select(-source) %>% distinct()
write_csv(bigUniqueFindTest2, "intraCompare2_0-2-15_FullData.csv")

head(bigUniqueFind %>% 
  dplyr::select(gene, cluster1, cluster2) %>%
  filter(
    cluster1 %in% c(0, 2, 15),
    !(cluster2 %in% c(0, 2, 15))) %>% 
  count(gene, cluster1) %>% filter(n == max(n)))

bigUniqueFindTest2 %>% distinct() %>% group_by(gene, cluster1, cluster2) %>% count()

# Used HOMER to find predicted FXR binding motifs----
FXR_sites <- read_tsv("FXR_regions.bed",
  col_names = c("Chromosome", "Start", "End", "Name", "Score", "Strand", "Chrom2", "Start2", "End2",
    "ensembl_transcript_id_version", "Distance"))

FXR_sites <- FXR_sites %>% dplyr::select(10, 11) %>% left_join(bm_list)
FXR_sites <- FXR_sites %>% dplyr::select(2,3,4)
FXR_sites <- distinct(FXR_sites)

FXR_5kb_tally <- FXR_sites %>%
  filter(Distance < 10000) %>%
  mutate(avgDist = mean(Distance)) %>% add_count(external_gene_name)


FXR_sites <- read_tsv("FXR_regions.bed",
  col_names = c("Chromosome", "Start", "End", "Name", "Score", "Strand", "Chrom2", "Start2", "End2",
    "ensembl_transcript_id_version", "Distance"))

FXR_sites <- FXR_sites %>% dplyr::select(10, 11) %>% left_join(bm_list)
FXR_sites <- FXR_sites %>% dplyr::select(2,3,4)
FXR_sites <- distinct(FXR_sites)

FXR_5kb_tally <- FXR_sites %>%
  filter(Distance < 10000) %>%
  mutate(avgDist = mean(Distance)) %>% add_count(external_gene_name)


# Making gene lists for GO Term Functional Analysis with DAVID or Metascape----
# Make the "Foreground" Gene lists
bm_list <- bm_list %>% group_by(ensembl_gene_id, external_gene_name) %>%
  nest()

for (i in 0:26) {
  FAmarkers %>% filter(cluster == i & (pct.1 - pct.2) > 0.05 & avg_logFC > log10(2)) %>%
    arrange(p_val_adj, desc(avg_logFC)) %>%
    left_join(bm_list, by = c("gene" = "external_gene_name")) %>% na.omit() %>%
    dplyr::select("ensembl_gene_id") %>%
    write_tsv(paste0("Output/FindAllMarkers/FindAllMarkers_Cluster_", i, ".tsv"),
      col_names = FALSE)
}
df <- tibble()

for (i in 0:26) {
  temp <- FAmarkers %>% filter(cluster == i & (pct.1 - pct.2) > 0.05 & avg_logFC > log10(2)) %>%
    arrange(p_val_adj, desc(avg_logFC)) %>%
    left_join(bm_list, by = c("gene" = "external_gene_name")) %>% na.omit() %>%
    dplyr::select(i = "ensembl_gene_id")
  df[i + 1] <- temp
}
write_tsv(df, "Output/FindAllMarkers/FindAllMarkers_AllClusters_colbind.tsv")



FAmarkers %>% filter((pct.1 - pct.2) > 0.05, avg_logFC > log10(2)) %>%
  group_by(cluster) %>% count()

bigUniqueFind %>% 
  filter(cluster2 == 0 & pct.1 > 0 & avg_logFC > log10(2) & !(cluster1 %in% c(20, 25, 26))) %>%
  top_n(1000, wt = avg_logFC) %>% dplyr::select(gene) %>% 
    left_join(bm_list, by = c("gene" = "external_gene_name"))

for (i in 0:26L) {
  bigUniqueFind %>%
    filter(cluster2 == i & pct.1 > 0 & avg_logFC > log10(2) & !(cluster1 %in% c(20, 25, 26))) %>%
    left_join(bm_list, by = c("gene" = "external_gene_name")) %>% na.omit() %>%
    top_n(1000, wt = avg_logFC) %>%
    dplyr::select("ensembl_gene_id") %>%
    write_tsv(paste0("Output/FindAllMarkers/FindMarkers_Cluster_", i, "_Background_1K.tsv"),
      col_names = FALSE)
}

for (i in 0:26L) {
  print(paste("Filtering Cluster", i, "... number of genes:", sep = " "))
  
  rdTest <- refMarkers %>%
    filter(pct.1 > pct.2 & avg_logFC > log10(5) & cluster == i & !(gene %in% riboEnsembl$external_gene_name)) %>%
    arrange(p_val_adj, desc(avg_logFC)) %>%
    left_join(bm_list, by = c("gene" = "external_gene_name")) %>% na.omit() %>%
    dplyr::select("ensembl_gene_id", "gene", "cluster") %>%
    left_join(dplyr::select(bigUniqueFind, gene, cluster2), by = "gene")
  
  selGenes <- filter(rdTest, cluster2 == cluster)$gene %>% unique()
  
  rdTest <- rdTest %>%
    filter(!(gene %in% selGenes)) %>%
    dplyr::select(2,3) %>%
    distinct()
  
  print(nrow(rdTest))
  
  rdTest %>%
    dplyr::select(gene) %>%
    write_tsv(paste0("Output/FindAllMarkers/FindMarkers_Cluster_", i, "_UniqueFiltered.noRPGs.tsv"),
      col_names = FALSE)
}

Cluster5_6 <- FindMarkers(df_integrated, ident.1 = c(5,6))

rdTest <- Cluster5_6 %>%
  rownames_to_column(var = "gene") %>%
    filter(pct.1 > pct.2 & avg_logFC > log10(5) & !(gene %in% riboEnsembl$external_gene_name)) %>%
    arrange(p_val_adj, desc(avg_logFC)) %>%
    left_join(bm_list, by = c("gene" = "external_gene_name")) %>% na.omit() %>%
    dplyr::select("ensembl_gene_id", "gene") %>%
    left_join(dplyr::select(bigUniqueFind, gene, cluster2), by = "gene")
  
  selGenes <- filter(rdTest, cluster2 %in% c(5,6))$gene %>% unique()
  
  rdTest <- rdTest %>%
    filter(!(gene %in% selGenes)) %>%
    dplyr::select(2) %>%
    distinct()

geneList <- dfExpression$gene[str_which(dfExpression$gene,"^(rpl|rps)")]
for (i in 1:length(geneList)) {
  FeaturePlot(df_integrated, features = geneList[i], label = TRUE)
  ggsave(paste0("Output/Plots/FeaturePlots/APOA/FeaturePlot_", geneList[[i]], ".png"),
    height = 5, width = 5)
}

Cluster1_7_vs10 <- FindMarkers(df_integrated, ident.1 = c(1, 7), ident.2 = 10)

rdTest <- Cluster1_7_vs10 %>%
  rownames_to_column(var = "gene") %>%
    filter(pct.1 > pct.2 & avg_logFC > log10(5) & !(gene %in% riboEnsembl$external_gene_name)) %>%
    arrange(p_val_adj, desc(avg_logFC)) %>%
    left_join(bm_list, by = c("gene" = "external_gene_name")) %>% na.omit() %>%
    dplyr::select("ensembl_gene_id", "gene") %>%
    left_join(dplyr::select(bigUniqueFind, gene, cluster2), by = "gene")

selGenes <- filter(rdTest, cluster2 %in% c(5,6))$gene %>% unique()
  
rdTest <- rdTest %>%
  filter(!(gene %in% selGenes)) %>%
  dplyr::select(2) %>%
  distinct()


Clusters4_16 <- FindMarkers(object = df_integrated,
                       ident.1 = "MUT",
                       group.by = "library",
                       subset.ident = c(4, 16)) %>%
  rownames_to_column(var = "gene")

write_csv(Clusters4_16, "Output/DiffExpTests/MUTvsWT/WilcoxRankSum/MutvWT_Cluster_4-16_FinalUMAP.csv")
