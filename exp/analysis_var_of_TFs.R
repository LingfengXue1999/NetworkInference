# variability of ATAC-seq data  # 2021/8/16
# load data
library(openxlsx)
data_list = list()
cells = c("H1ESC","K562","GM12878","EML","mESC","TF1","HL60","BJ")
species = c(rep("human",3),rep("mouse",2),rep("human",3))
for(i in 1:length(cells)){
  a = read.xlsx("./ATACseq_var/2015_accessibility_variability_data.xlsx",sheet = i)
  data_list[[i]] = a
}
# transform list to dataframe
Anno = data_list[[1]]$Annotation
for(i in 1:length(data_list)){Anno = intersect(Anno,data_list[[i]]$Annotation)}
mat = matrix(data = NA, nrow = length(Anno), ncol = length(cells))
for(i in 1:length(data_list)){
  tmp = data_list[[i]]
  mat[,i] = tmp$Variability[match(Anno,tmp$Annotation)]
}
rownames(mat) = Anno
colnames(mat) = cells
# extract genes
genes = Anno[grep("motif",Anno)]
mat = mat[genes,]
# change names
# detect "-"
str_times= sapply(genes,function(x){return(length(strsplit(x,"-")[[1]])-1)})
TF_list = list()
for(i in 1:length(genes)){
  TFs = strsplit(genes[i],"-")[[1]][1:str_times[i]]
  TF_list[[length(TF_list)+1]] = TFs
}
TF_list[[12]] = "NANOG"
TF_list[[13]] = "OCT4"
TF_list[[14]] = "SMAD1"
# calculate mean variability for each motif
var_m = apply(mat,1,mean)
# calcualte all TFs
TFs_all = NULL
for(i in 1:length(TF_list)){
  TFs_all = c(TFs_all,TF_list[[i]])
}
TFs_all = unique(TFs_all)
# map motif to TFs
TF_var = list()
for(i in 1:length(TFs_all)){
  TF_var[[i]] = NA
  names(TF_var)[i] = TFs_all[i]
}
for(i in 1:length(var_m)){
  TFs = TF_list[[i]]
  for(TF in TFs){
    if(is.na(TF_var[TF])){
      TF_var[[TF]]= var_m[i]
    }
    else{
      TF_var[[TF]] = c(TF_var[[TF]],var_m[i])
    }
  }
}
# calculate TF variability for all motifs
var_m_TF = NULL
for(i in 1:length(TF_var)){
  var_m_TF = c(var_m_TF,mean(TF_var[[i]]))
}
names(var_m_TF) = TFs_all
names(var_m_TF) = toupper(TFs_all)
# save results
# save(var_m_TF,file = "./ATACseq_var/processed_vector.RData")



tmp = as.character(sapply(genes,function(x){return(strsplit(x,"-")[[1]][1])}))
tmp[12:14] = c("NANOG","OCT4","SMAD1")
rownames(mat) = genes = toupper(tmp)

# visualize data -- pca
if(T){
  a = t(mat)    # genes as variables
  pr.out <- prcomp (a,scale. = T)
  library(ggplot2)
  project <- predict(pr.out, newdata=a)[,1:2]
  project.new <- cbind(as.data.frame(project), cluster=as.factor(cells))
  head(project.new)
  ggplot(project.new, aes(x=PC1, y=PC2,color = cluster)) + geom_point() +geom_text(aes(label=cluster), hjust=0.5, vjust=1)
}
var_mat_mouse = mat[,which(species=="mouse")]
var_mat_human = mat[,which(species=="human")]

# save results
# var_mat = mat
# save(var_mat,var_mat_mouse,var_mat_human,file = "./ATACseq_var/processed_matrix.RData")

# distribution
boxplot(mat)
apply(mat,2,function(x){return(quantile(x,0.8))})

thr = 0.1
var_TFs = NULL
for(i in 1:ncol(mat)){
  tmp = mat[,i]
  var_TFs[[i]] = rownames(mat)[which(tmp>quantile(tmp,1-thr))]
}
var_TFs_sum = unlist(var_TFs)
sort(table(var_TFs_sum),decreasing = T)
# plot correlation between results
cor(mat)
# correlation between hFB results and average
# load hFB data
load(file = "./ATACseq_var/forebrain/processed_var_vec.RData")
var.hFB = var_vec
overlap_genes = intersect(names(var.hFB),rownames(mat))
mat.sub = mat[overlap_genes,]
mat.com = cbind(mat.sub,var.hFB[overlap_genes])
colnames(mat.com) = c(colnames(mat),"hFB")
cor(mat.com)
