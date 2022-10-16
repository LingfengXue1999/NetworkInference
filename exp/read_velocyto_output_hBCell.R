## read velocyto output 2020/8/18
#density
density.matrix = function(a){
  a = as.matrix(a)
  t = as.logical(a[,])
  sum = sum(t)
  return(sum/length(a))
}
density.sparsematrix = function(a){
  nonzero = length(a@x)
  all = dim(a)[1]*dim(a)[2]
  return(nonzero/all)
}



dir = './more_datasets_0301/hBCell'
library(Matrix)
# preprocess velocyto output
if(T){
  
  s = readMM(paste(dir,'/spliced_matrix.mtx',sep = ''))
  u = readMM(paste(dir,'/unspliced_matrix.mtx',sep = ''))
  
  genes = read.table(paste(dir,'/genes.txt',sep = ''),sep = '\t')
  cells = read.table(paste(dir,'/cells.txt',sep = ''),sep = '\t')
  
  genes = as.character(as.matrix(genes))
  cells = as.character(as.matrix(cells))

  rownames(u) = rownames(s) = genes
  colnames(u) = colnames(s) = cells
  
  ##filter total-zero genes
  rsum = rowSums(u) + rowSums(s)
  index = rsum>0
  sum(index)
  u.mat = u[index,]
  s.mat = s[index,]
  genes.sub = genes[index]
  rm(list = c('u','s','genes'))
  
  
  
  density.sparsematrix(u.mat)     ##4.1%
  density.sparsematrix(s.mat)     ##8.9%
  
  save(u.mat,s.mat, file = paste(dir,'/unspliced_spliced_matrix.RData',sep = ''))
  
}

load(file =  paste(dir,'/unspliced_spliced_matrix.RData',sep = ''))

if(T){
  genes.mat = rownames(u.mat)
  #intron names
  intron_names = paste(genes.mat,'-intron',sep = '')
  intron_names[1:10]
  rownames(u.mat) = intron_names
  
  #combind intron data with exon data
  data_all_mat = rbind(s.mat,u.mat)
  
}

#Normalize
if(T){
  library(Seurat)
  data.seurat <- CreateSeuratObject(counts = data_all_mat)  
  dim(data.seurat)  ##  
  data.seurat <- NormalizeData(data.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  
  
  ####QC
  #remove low expression genes: only for total counts (exprMat)
  data.seurat.1 = CreateSeuratObject(data.seurat[["RNA"]]@data, min.cells = 5)     ####parameter: at least expressed in 5 cells

  dim(data.seurat.1)   ## 33388 13301
  rm('data.seurat')
  rm('data_all_mat')
  #subsetMT genes
  VlnPlot(data.seurat.1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)    
  FeatureScatter(data.seurat.1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
   
  # select cells: nCount_RNA>3000
  data.sub = subset(data.seurat.1,subset = nCount_RNA > 1500)
  cell.names.sub = colnames(data.sub[["RNA"]]@data)     #selected cells      ## 37767 18994
  dim(data.sub)
  VlnPlot(data.sub, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)    
  rm('data.seurat.1')
  
  # save filtered results
  save(u.mat,s.mat,cell.names.sub,file = paste(dir,'/unspliced_spliced_matrix.RData',sep = ''))
  
  a = data.sub[["RNA"]]@data
  genes.exon = genes.mat
  genes.intron = intron_names
  intron = a[grepl('-intron',rownames(a)),] 
  exon = a[!grepl('-intron',rownames(a)),]
  genes.intron = sub('-intron','', rownames(intron))   # remove '-intron'
  rownames(intron) = genes.intron
  # intron = as.matrix(intron)          # 14103  1720
  # exon = as.matrix(exon)              # 13558  1720
  dim(intron)
  dim(exon)
  density.sparsematrix(intron)              # 6.7%   
  density.sparsematrix(exon)                # 14%  
}

# investigate the dataset
if(T){
  # highly variable genes
  data.sub = FindVariableFeatures(data.sub, selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(data.sub)
  data.sub <- ScaleData(data.sub, features = VariableFeatures(object = data.sub))
  # PCA
  data.sub <- RunPCA(data.sub, features = VariableFeatures(object = data.sub))
  #print(data.sub[["pca"]], dims = 1:5, nfeatures = 5)		
  #VizDimLoadings(data.sub, dims = 1:2, reduction = "pca")
  DimPlot(data.sub, reduction = "pca")		
  
  # UMAP
  data.sub <- RunUMAP(data.sub, dims = 1:10)
  DimPlot(data.sub, reduction = "umap" )
  
  
}

# clustering
if(T){
  data.sub <- FindNeighbors(data.sub, dims = 1:10)
  data.sub <- FindClusters(data.sub, resolution = 0.2)
  DimPlot(data.sub, reduction = "umap")
}


##further QC
thr = 0.1
# genes = which(rowSums(u.mat[,cell.names.sub])>thr*length(cell.names.sub))
# intron.sub = intron[sub('-intron','', names(genes)) ,]
# density.sparsematrix(intron.sub)  # 21%
a = intron
t = matrix(as.logical(a),nrow = nrow(a))
genes = rowSums(t)

index = genes[]>thr*ncol(intron)    # 3504 intron
sum(index)
intron.sub = intron[index,]
density.matrix(intron.sub)     # 23%

# genes = which(rowSums(s.mat[,cell.names.sub])>thr*length(cell.names.sub))
# exon.sub = exon[intersect(names(genes),rownames(exon)),]
# density.sparsematrix(exon.sub)   # 39%
a = exon
t = matrix(as.logical(a),nrow = nrow(a))
genes = rowSums(t)

index = genes[]>thr*ncol(intron)     # 6293 exons remained
sum(index)
exon.sub = exon[index,]
density.matrix(exon.sub)     # 33%

## get TF list from RcisTarget package: ## mouse
library(RcisTarget)
data(motifAnnotations_hgnc)
TF_list = unique(motifAnnotations_hgnc$TF)

TFs = intersect(rownames(exon.sub),TF_list)
length(TFs)    # 612 TFs
exon.sub = as.matrix(exon.sub)
intron.sub = as.matrix(intron.sub)

save(exon.sub, intron.sub, TFs, file = paste(dir,'/intron_exon_data_processed.RData',sep = ''))




