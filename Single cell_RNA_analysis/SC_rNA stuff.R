

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(SCINA)
max.overlaps=Inf







# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/filtered_gene_bc_matrices/IHNV_1d_sc/")
DATA <- pbmc.data
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "IHNV", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
PercentageFeatureSet(pbmc, pattern = "^mt-")
head(pbmc@meta.data, 18)
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc[["RNA"]]@scale.data
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:19)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 30)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:30)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
# find all markers of cluster 2
FindNeighbors(pbmc, 1:19)
cluster2.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 12), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
load(system.file('filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/barcodes.tsv','filtered_gene_bc_matrices/filtered_gene_bc_matrices/IHNV_1d_sc/', package = "SCINA"))
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = c("Immune system","Brain") # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 
pbmc@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pbmc@meta.data$customclassif[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("zgc:153704", "zgc:153704",   "cd59", "ccl34a.4" , "mbpa",         
                           "BX649485.4",     "aplnra",        "CABZ01101977.1", "id1",           
                            "plp1b","tuba8l3" ))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("PACAP", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
VlnPlot(pbmc,features="nFeature_RNA")
FeaturePlot(pbmc, features = c("zgc:153704", "zgc:153704",   "cd59", "ccl34a.4" , "mbpa",         
                               "BX649485.4",     "aplnra",        "CABZ01101977.1", "id1",           
                               "plp1b","tuba8l3" ))
      


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
FeaturePlot(pbmc, features = c("ch211-214p16.1",
                                  "ccl38.6",          
                                     "lgals3bpb",       
                                     "si:busm1-266f07.2",
                                     "lgals9l1",        
                                     "lta4h",           
                                    "cort",           
                                     "nrgna",           
                                   "CABZ01101977.1",   
                                      "pyyb"  ))
FeaturePlot(pbmc, features = c( "ifih1" , "tlr3" , "irf7" , "tbk1" , "stat1" ,"mx1", "isg15" , "rsad2"))

top30 <- head(VariableFeatures(pbmc), 30)
top30
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
plot2
#this is just making the feature plots of the top 30 genes 
FeaturePlot(pbmc, features = c("zgc:153704", "cd59", "ccl34a.4", "mbpa", "BX649485.4", 
                               "aplnra", "CABZ01101977.1", "id1", "plp1b", "tuba8l3", 
                               "g0s2", "olig1", "adcyap1b", "zgc:165461", "si:ch211-286c4.6", 
                               "BX649485.2", "olig2", "mpz", "wu:fj16a03", "pyyb", 
                               "ighv1-4", "BX005103.1", "flj13639", "igfbp1a", "epd", 
                               "fabp7a", "BX908782.2", "fabp11a", "gfap", "zgc:112332"))
FeaturePlot(pbmc, features = c("zgc:153704", "cd59", "ccl34a.4", "mbpa", "BX649485.4", 
                               "aplnra", "CABZ01101977.1", "id1", "plp1b", "tuba8l3"))
FeaturePlot(pbmc, features = c("g0s2", "olig1", "adcyap1b", "zgc:165461", "si:ch211-286c4.6", 
                               "BX649485.2", "olig2", "mpz", "wu:fj16a03", "pyyb"))
FeaturePlot(pbmc, features = c("ighv1-4", "BX005103.1", "flj13639", "igfbp1a", "epd", 
                               "fabp7a", "BX908782.2", "fabp11a", "gfap", "zgc:112332"))

