library(Seurat)

files=list.files('./GSE135927_RAW/','^GSM')
print(files)
samples=strsplit(files,'_')

#extract the datasets into two folders, then rename to barcodes, 
#features, and matrix for subsequent seurat object building
lapply(unique(samples),function(x){
  filename = paste(x[1], x[2], x[3], sep = "_")
  foldername = paste("GSE135927_RAW/", x[1], sep = "")
  if (!dir.exists(foldername)){
    dir.create(foldername,recursive = T)
  }
  file.rename(file.path("GSE135927_RAW", filename), file.path(
    "GSE135927_RAW", x[1], x[3]))
})

#create seurat objects for both samples
samples=list.files("GSE135927_RAW/")
seuratList = lapply(samples, function(y){
  folderpath = file.path("GSE135927_RAW", y)
  CreateSeuratObject(counts = Read10X(folderpath), project = y)
})

#merge the seruatList into a Large Seurat object
#https://satijalab.org/seurat/articles/merge_vignette.html
seurat.combined <- merge(seuratList[[1]],y = c(seuratList[[2]]), add.cell.ids = samples,
  project = "Large_Seurat_Object")

table(seurat.combined$orig.ident)#track which sample originated from

##checkpt
save(seurat.combined,file = 'seurat.combined.Rdata')
seurat.combined=get(load(file = 'seurat.combined.Rdata'))

#Quality control and making of violin plots
#reference from https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html
seurat.combined[["percent.mt"]] <-
  PercentageFeatureSet(seurat.combined, pattern = "^MT-")
seurat.combined[["percent.ribo"]] <-
  PercentageFeatureSet(seurat.combined, pattern = "^RP")
plot1 <-
  FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <-
  FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(
  seurat.combined,
  features = c("percent.ribo", "percent.mt"),
  ncol = 2
)
VlnPlot(
  seurat.combined,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2
)
VlnPlot(
  seurat.combined,
  features = c("percent.ribo", "nCount_RNA"),
  ncol = 2
)
#TODO: interpret violin plots to generate filtering parameters

# Each sample, from a total of 1359 naive and 683 SR1 cells, was
# first filtered to remove cells with low gene counts (naive, n = 92;
# SR1, n = 37) that arise from aborted sequencing, and gene
# expression was normalized between cells.
##due to unknown parameter for filtering, parameters used are based on
##https://blog.katastros.com/a?ID=af0c80dc-5f37-4a17-b7b4-31edb9e51c2b
filtered.seurat.combined <-
  subset(seurat.combined, subset = nFeature_RNA > 200 &
           nCount_RNA > 1000 &
           percent.mt < 20)
filtered.seurat.combined
table(filtered.seurat.combined$orig.ident)

##checkpt
save(filtered.seurat.combined, file = 'filtered.seurat.combined.Rdata')
filtered.seurat.combined = get(load(file = 'filtered.seurat.combined.Rdata'))

#normalize data
filtered.seurat.combined <-
  NormalizeData(
    filtered.seurat.combined,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

# Afterward, variable expression of genes was determined.
# Naive and SR1 samples were then merged and aligned
# dimensions for t - distributed stochastic neighbor
# embedding were calculated to identify unique cell clusters,
# and the differential expression of genes between cell clusters was determined.
##workflow follows https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html
filtered.seurat.combined <-
  FindVariableFeatures(filtered.seurat.combined,
                       selection.method = "vst",
                       nfeatures = 2000)

top10 <- head(VariableFeatures(filtered.seurat.combined), 10)
#feature scatter plot 
plot1 <- VariableFeaturePlot(filtered.seurat.combined)
plot2 <- LabelPoints(plot = plot1,
                     points = top10,
                     repel = FALSE)
plot1
plot2

#scaling
all.genes <- rownames(filtered.seurat.combined)
filtered.seurat.combined <-
  ScaleData(filtered.seurat.combined, features = all.genes)

#PCA analysis
filtered.seurat.combined <-
  RunPCA(filtered.seurat.combined,
         features = VariableFeatures(object = filtered.seurat.combined))
print(filtered.seurat.combined[["pca"]],
      dims = 1:5,
      nfeatures = 5)
VizDimLoadings(filtered.seurat.combined,
               dims = 1:2,
               reduction = "pca")
DimPlot(filtered.seurat.combined, reduction = "pca")
DimHeatmap(
  filtered.seurat.combined,
  dims = 1:5,
  cells = 500,
  balanced = TRUE
)

#determine dimensionality
ElbowPlot(filtered.seurat.combined)#elbow at 18

#parameters from https://blog.katastros.com/a?ID=af0c80dc-5f37-4a17-b7b4-31edb9e51c2b
filtered.seurat.combined <-
  FindNeighbors(filtered.seurat.combined, dims = 1:18)
filtered.seurat.combined <-
  FindClusters(filtered.seurat.combined, resolution = 0.9)

table(filtered.seurat.combined@meta.data$RNA_snn_res.0.9)


#t-SNE clustering visualisation
set.seed(123)
filtered.seurat.combined <-
  RunTSNE(object = filtered.seurat.combined,
          dims = 1:18,
          do.fast = TRUE)
DimPlot(filtered.seurat.combined,
        reduction = "tsne",
        label = T)
DimPlot(
  filtered.seurat.combined,
  reduction = "tsne",
  label = T,
  split.by = 'orig.ident'
)#seperate by samples

##checkpoint
saveRDS(filtered.seurat.combined, file = "./progenitorTcells.rds")

#further analysis: clusters cell type labelling
seurat.markers <-
  FindAllMarkers(filtered.seurat.combined,
                 min.pct = 0.25,
                 min.diff.pct = 0.25)
ref <- MonacoImmuneData()#bulk RNA-seq data of sorted human immune cells
seurat_for_SingleR <- GetAssayData(filtered.seurat.combined, slot = "data")
clusters <- filtered.seurat.combined@meta.data$seurat_clusters
#use previously determined cluster and annotate them
celltype_pred <-
  SingleR(
    seurat_for_SingleR,
    ref,
    labels = ref$label.fine, #fine for immune cells 
    clusters = clusters,
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts"
  )
celltype = data.frame(
  ClusterID = rownames(celltype_pred),
  celltype = celltype_pred$labels,
  stringsAsFactors = F
)

#add cell type to seruat object 
for (i in 1:nrow(celltype)) {
  filtered.seurat.combined@meta.data[which(
    filtered.seurat.combined@meta.data$seurat_clusters == 
      celltype$ClusterID[i]), 'celltype'] <-
    celltype$celltype[i]
}
plot = DimPlot(
  filtered.seurat.combined,
  group.by = "celltype",
  label = T,
  label.size = 5,
  reduction = 'tsne'
)
plot

#checkpoint
saveRDS(filtered.seurat.combined, file = "./progenitorTcells_final.rds")


