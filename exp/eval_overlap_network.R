# evaluate overlap network # 2021/9/9
# load NI result
sample_ids = c(1,15)
sample_id = 1
mat.list = list()
for(sample_id in sample_ids){
  # load network
  if(T){
    sample_names = c("hForebrain","BoneMarrow","epithelium","DentateGyrus","endocrinogenesis","hESC_0,2",
                     "hESC_day0","hESC_day2","hESC_day6_Mesen","hESC_day6_Endo","hESC_day4","hESC_day8_Mesen",
                     "hESC_2018","hESC_2019","mND")
    sample_sources = c("human","mouse","mouse","mouse","mouse",rep("human",9),"mouse")
    sample_name = sample_names[sample_id]
    sample_source = sample_sources[sample_id]
    if(sample_id==1){
      # load exon result
      load(file = "./10X_datasets/unspliced_spliced_matrix/Velocyto/Forebrain/network_inference_result/genie3_exon_Velo_hFB.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./10X_datasets/unspliced_spliced_matrix/Velocyto/Forebrain/network_inference_result/genie3_inex_Velo_hFB.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==2){
      # load exon result
      load(file = "./10X_datasets/unspliced_spliced_matrix/Velocyto/BoneMarrow/network_inference_result/genie3_exon_Velo_BM.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./10X_datasets/unspliced_spliced_matrix/Velocyto/BoneMarrow/network_inference_result/genie3_inex_Velo_BM.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==3){
      # load exon result
      load(file = "./10X_datasets/unspliced_spliced_matrix/Velocyto/epithelium/network_inference_result/genie3_exon_Velo_Ep.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./10X_datasets/unspliced_spliced_matrix/Velocyto/epithelium/network_inference_result/genie3_inex_Velo_Ep.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==4){
      # load exon result
      load(file = "./10X_datasets/unspliced_spliced_matrix/scVelo/DentateGyrus/network_inference_result/genie3_exon_all2.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./10X_datasets/unspliced_spliced_matrix/scVelo/DentateGyrus/network_inference_result/genie3_inex_all2.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==5){
      # load exon result
      load(file = "./10X_datasets/unspliced_spliced_matrix/scVelo/endocrinogenesis/network_inference_result/genie3_exon_6.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./10X_datasets/unspliced_spliced_matrix/scVelo/endocrinogenesis/network_inference_result/genie3_inex_6.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==6){
      # load exon result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_exon_hESC_day0_and_day2_1000cells.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_inex_hESC_day0_and_day2_1000cells.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==7){
      # load exon result
      load(file = "./RNAseq/hESC/SRR9117953_2019/genie3_exon_output_2019.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/SRR9117953_2019/genie3_inex_output_2019.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==8){
      # load exon result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_exon_hESC_day2_1500cells.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_inex_hESC_day2_1500cells.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==9){
      # load exon result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_exon_hESC_day6_Mesenchymal.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_inex_hESC_day6_Mesenchymal.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==10){
      # load exon result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_exon_hESC_day6_Endothelial.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_inex_hESC_day6_Endothelial.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==11){
      # load exon result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_exon_hESC_day4.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_inex_hESC_day4.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==12){
      # load exon result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_exon_hESC_day8_Mesenchymal.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/GSE131736/network_inference_result/genie3_inex_hESC_day8_Mesenchymal.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==13){
      # load exon result
      load(file = "./RNAseq/hESC/SRR6328624_2018/genie3_exon_output.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/SRR6328624_2018/genie3_inex_output.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==14){
      # load exon result
      load(file = "./RNAseq/hESC/SRR9117953_2019/genie3_exon_output_2019.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/hESC/SRR9117953_2019/genie3_inex_output_2019.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    if(sample_id==15){
      # load exon result
      load(file = "./RNAseq/GSE152125/genie3_exon_mND.RData")
      exon.mat = weightMatrix
      
      # load inex result
      load(file = "./RNAseq/GSE152125/genie3_inex_mND.RData")
      mat = weightMatrix
      # extract intron targets
      targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
      mat = mat[,targets]
      # change intron names
      colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
      # remove self regulation
      for(i in 1:nrow(mat)){
        if(is.element(rownames(mat)[i],colnames(mat))){
          mat[i,rownames(mat)[i]]=0
        }
      }
      inex.mat = mat
    }
    # check NaN
    nan.exon = which(is.nan(exon.mat[1,]))
    nan.inex = which(is.nan(inex.mat[1,]))
    if(length(nan.exon)!=0){
      exon.mat = exon.mat[,-nan.exon]
    }
    if(length(nan.inex)!=0){
      nan.inex = nan.inex[,-nan.inex]
    }
    dim(exon.mat)
    dim(inex.mat)
  }
  # select overlap genes
  overlap_targets = intersect(colnames(exon.mat),colnames(inex.mat))
  exon.mat = exon.mat[,overlap_targets]
  inex.mat = inex.mat[,overlap_targets]
  dim(exon.mat)
  rownames(exon.mat) = rownames(inex.mat) = toupper(rownames(exon.mat))
  colnames(exon.mat) = colnames(inex.mat) = toupper(colnames(exon.mat))
  
  mat.list[[length(mat.list)+1]] = list(exon = exon.mat,inex = inex.mat)
}
TFs = list(hFB = rownames(mat.list[[1]][["exon"]]),mFB = rownames(mat.list[[2]][["exon"]]))
overlap_TFs = intersect(TFs[[1]],TFs[[2]])
targets = list(hFB = colnames(mat.list[[1]][["exon"]]),mFB = colnames(mat.list[[2]][["exon"]]))
overlap_targets = intersect(targets[[1]],targets[[2]])

hFB.in = mat.list[[1]][["exon"]][overlap_TFs,overlap_targets]
hFB.ex = mat.list[[1]][["inex"]][overlap_TFs,overlap_targets]
mFB.in = mat.list[[2]][["exon"]][overlap_TFs,overlap_targets]
mFB.ex = mat.list[[2]][["inex"]][overlap_TFs,overlap_targets]

# compare overlap between two datasets
cal_overlap = function(vec1,vec2,n_top){
  targets1 = names(sort(vec1,decreasing = T))[1:n_top]
  targets2 = names(sort(vec2,decreasing = T))[1:n_top]
  overlap_targets = intersect(targets1,targets2)
  return(overlap_targets)
}
n_top = 500
df.eval = data.frame()
for(i in 1:length(overlap_TFs)){
  TF = overlap_TFs[i]
  vec1 = hFB.in[TF,]
  vec2 = mFB.in[TF,]
  over_targets = cal_overlap(vec1,vec2,n_top)
  n1 = length(over_targets)
  p1 = n1/n_top
  
  vec1 = hFB.ex[TF,]
  vec2 = mFB.ex[TF,]
  over_targets = cal_overlap(vec1,vec2,n_top)
  n2 = length(over_targets)
  p2 = n2/n_top
  
  
  df.eval = rbind(df.eval,c(i,n_top,n1,p1,n2,p2))
}
colnames(df.eval) = c("TF_id","n_top","n_over_in","ratio_in","n_over_ex","ratio_ex")
rownames(df.eval) = overlap_TFs
random_ratio = n_top/length(overlap_targets)

pdf(file = paste("hFB_mFB_consistency.pdf",sep = ""), width = 6,height = 6)

boxplot(list(a = df.eval$ratio_in,b = df.eval$ratio_ex), names = c("intron","exon"),ylab = "overlap ratio",
        main = "Compare consistency between hFB and mFB")
abline(h = random_ratio, col = "grey")
wilcox.test(df.eval$ratio_in,df.eval$ratio_ex,alternative = "greater")$p.value
dev.off()


boxplot(df.eval$ratio_in/df.eval$ratio_ex,ylab = "overlap ratio: intron/exon",
        main = "compare intron with exon for each TF")
abline(h = 1, col = "grey")

plot(df.eval$ratio_ex,df.eval$ratio_in)
abline(a = 0, b = 1, col = "green")

genesets = c("DLX1","DLX2","GSX1","GSX2","ASCL1","OLIG2","ARX","GAD1","MAF",
             "DLX5","DLX6","PBX1","ETV1","ER81","SIX3","ZFHX1B","FOXP4",
             "SP9","SP8","VAX1","LHX8","LHX6","CXCR7","SOX6","MAFB","GBX1",
             "GBX2","SHH","HELT","SIX3","EBF1","HES5","OTP")
over_genes = intersect(overlap_TFs,genesets)
View(df.eval[over_genes,])


# enrichment analysis
TF = "NEUROD2"
vec1 = hFB.in[TF,]
vec2 = mFB.in[TF,]
over_targets1 = cal_overlap(vec1,vec2,n_top)

vec1 = hFB.ex[TF,]
vec2 = mFB.ex[TF,]
over_targets2 = cal_overlap(vec1,vec2,n_top)

over_targets1
over_targets2
library(AnnotationDbi)
library(org.Hs.eg.db) 
geneset1 = mapIds(org.Hs.eg.db,over_targets1,keytype = "SYMBOL",column = "ENTREZID")
geneset2 = mapIds(org.Hs.eg.db,over_targets2,keytype = "SYMBOL",column = "ENTREZID")
# GO analysis   -- BP
ego1 <- enrichGO(OrgDb="org.Hs.eg.db", gene = geneset1, ont = "BP", pvalueCutoff = 0.05, readable= TRUE)
ego2 <- enrichGO(OrgDb="org.Hs.eg.db", gene = geneset2, ont = "BP", pvalueCutoff = 0.05, readable= TRUE)
dotplot(ego,showCategory=5,title="Enrichment GO Top5") 
library("cowplot")
pdf(file = "./GO_NEUROD2_overlap_targets.pdf",width = 16,height = 9)
p1 = barplot(ego1, showCategory=5,title="NEUROD2_overlap_targets:intron (n = 34)") 
p2 = barplot(ego2, showCategory=5,title="NEUROD2_overlap_targets:exon (n = 42)") 
plot_grid(p1,p2,nrow = 1,ncol = 2)
dev.off()


