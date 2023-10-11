# scType test
```
use this code:

library(knitr)
library(Seurat)
library(ggplot2)
library(HGNChelper)
library(dplyr)
options(Seurat.object.assay.version = "v4")
#source("C:/Users/admin/Desktop/sc-type-master/R/gene_sets_prepare.R")
#source("C:/Users/admin/Desktop/sc-type-master/R/sctype_score_.R")
source("./gene_sets_prepare.R")
source("./sctype_score_.R")
pbmc.data <- Read10X(data.dir = "C:/Users/admin/Desktop/filtered_feature_bc_matrix")
#pbmc.data <- Read10X(data.dir = "C:/Users/admin/Desktop/raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3,
    min.features = 200)
# normalize data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
# percent.mt < 5) # make some filtering based on QC metrics visualizations, see
# Seurat tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# cluster and visualize
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)


pbmc <- RunUMAP(pbmc, dims = 1:15)
#DimPlot(pbmc, reduction = "umap")

# DB file
db_ = "./ScTypeDB_fullM.xlsx"
tissue = "Liver"# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive,
                      gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq
# matrix.  In case Seurat is used, it is either pbmc[['RNA']]$scale.data
# (default), pbmc[['SCT']]$scale.data, in case sctransform is used for
# normalization, or pbmc[['integrated']]$scale.data, in case a joint analysis
# of multiple single-cell datasets is performed.
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

pbmc@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pbmc@meta.data$customclassif[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```
