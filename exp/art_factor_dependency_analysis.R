# Evaluate for each TF -- experimental datasets   2021/6/3
# evaluation function for all TFs
Evaluation_TF = function(mat,gt){
  df = NULL
  for(i in 1:nrow(mat)){
    weights = sort(mat[i,])
    gt_targets = which(gt[i,]==1)
    df.tmp = NULL
    for(j in 1:length(weights)){
      thr = weights[j]
      infer_targets = which(mat[i,]>=thr)
      true_targets = intersect(gt_targets,infer_targets)
      TP = length(true_targets)
      FP = length(infer_targets) - length(true_targets)
      FN = length(gt_targets) - length(true_targets)
      precision = TP/(TP+FP)
      recall = TP/(TP+FN)
      df.tmp = rbind(df.tmp,c(thr,precision,recall))
    }
    colnames(df.tmp) = c("thr","precision","recall")
    thr = df.tmp[,1]
    precision = df.tmp[,2]
    recall = df.tmp[,3]
    # calculat AUPR
    AUPR = 0
    for(k in 1:(length(thr)-1)){AUPR = AUPR+ (recall[k] - recall[k+1])*(precision[k]+precision[k+1])/2}
    random_precision = precision[1]
    # save data
    df = rbind(df,c(length(gt_targets),AUPR,random_precision))
  }
  rownames(df) = rownames(mat)
  colnames(df) = c("n_true_targets","AUPR","random_precision")
  df = as.data.frame(df)
  
  return(df)
}
Eval_top_for_TFs = function(mat,gt,top_num){
  df = NULL
  for(i in 1:nrow(mat)){
    weights = sort(mat[i,],decreasing = T)
    gt_targets = which(gt[i,]==1)
    top_precision = length(intersect(names(weights)[1:top_num],names(gt_targets)))/top_num
    random_precision = length(gt_targets)/length(weights)
    
    df = rbind(df,c(length(gt_targets),top_precision,random_precision,top_num))
  }
  rownames(df) = rownames(mat)
  colnames(df) = c("n_true_targets","top_precision","random_precision","top_num")
  df = as.data.frame(df)
  
  return(df)
}
list_to_mat = function(a){
  TFs = unique(a[,1])
  targets = unique(a[,2])
  mat = matrix(0,nrow = length(TFs),ncol = length(targets))
  rownames(mat) = sort(TFs)
  colnames(mat) = sort(targets)
  if(ncol(a)==2){
    for(i in 1:nrow(a)){
      mat[a[i,1],a[i,2]] = 1
    }
  }else{
    for(i in 1:nrow(a)){
      mat[a[i,1],a[i,2]] = as.numeric(a[i,3])
    }
  }
  return(mat)
}
# evaluation function for all targets
Evaluation_target = function(mat,gt){
  df = NULL
  for(i in 1:ncol(mat)){
    weights = sort(mat[,i])
    gt_TFs = which(gt[,i]==1)
    df.tmp = NULL
    for(j in 1:length(weights)){
      thr = weights[j]
      infer_TFs = which(mat[,i]>=thr)
      true_TFs = intersect(gt_TFs,infer_TFs)
      TP = length(true_TFs)
      FP = length(infer_TFs) - length(true_TFs)
      FN = length(gt_TFs) - length(true_TFs)
      TN = nrow(gt) - TP - FP - FN
      precision = TP/(TP+FP)
      recall = TP/(TP+FN)
      FPR = FP/(FP+TN)
      df.tmp = rbind(df.tmp,c(thr,precision,recall,FPR))
    }
    colnames(df.tmp) = c("thr","precision","recall","FPR")
    thr = df.tmp[,1]
    precision = df.tmp[,2]
    recall = df.tmp[,3]
    FPR = df.tmp[,4]
    # calculat AUPR
    AUPR = 0
    for(k in 1:(length(thr)-1)){AUPR = AUPR+ (recall[k] - recall[k+1])*(precision[k]+precision[k+1])/2}
    random_precision = precision[1]
    AUROC = 0
    for(k in 1:(length(thr)-1)){AUROC = AUROC+ (FPR[k] - FPR[k+1])*(recall[k]+recall[k+1])/2}
    
    # save data
    df = rbind(df,c(length(gt_TFs),AUPR,random_precision,AUROC))
  }
  rownames(df) = colnames(mat)
  colnames(df) = c("n_true_TFs","AUPR","random_precision","AUROC")
  df = as.data.frame(df)
  
  return(df)
}
Eval_top_for_targets = function(mat,gt,top_num){
  df = NULL
  for(i in 1:ncol(mat)){
    weights = sort(mat[,i],decreasing = T)
    gt_TFs = which(gt[,i]==1)
    top_precision = length(intersect(names(weights)[1:top_num],names(gt_TFs)))/top_num
    random_precision = length(gt_TFs)/length(weights)
    
    df = rbind(df,c(length(gt_TFs),top_precision,random_precision,top_num))
  }
  rownames(df) = colnames(mat)
  colnames(df) = c("n_true_TFs","top_precision","random_precision","top_num")
  df = as.data.frame(df)
  
  return(df)
}
library(ggplot2)

# load NI result 
sample_ids = c(1,15,2,5,3,4,13,14)
sample_ids = c(1,15,2,5)


sample_id = 1

# configurations for plot
# config1: effect of two parameters on targets, top precision, pvalue
## not evaluate TF; not evaluate proportion; not evaluate correlation between meancount and half-life; not show barplot
pdf("./effect_two_para_TP.pdf",width = 5,height = 5)   # plot for config1, top precision
pdf("./pvalue_for_4datasets.pdf",width = 7,height = 5) # plot for config1, pvalue
# config2: effect of two parameters on targets, barplot, only one sample
## par(mfrow = c(1,4)), plot barplot
pdf("./effect_two_para_barplot_hFB.pdf",width = 9,height = 4.5)
# config3: effect of mean count and sd on TFs
## not evaluate target; 
pdf("./effect_mean_TF.pdf",width = 5,height = 5)
# config4: effect of variability on TFs
## not evaluate TF; evaluate variability
pdf("./effect_var_TF_4datasets.pdf",width = 5,height = 5) # config4, 4 datasets
pdf("./effect_var_TF_summary.pdf",width = 5,height = 5)    # config4, all datasets summary
pdf("./effect_var_TF_hFBdata.pdf",width = 5,height = 5)   # config4, hFB as data
pdf("./effect_var_TF_summary_bar_pvalue.pdf",width = 5,height = 5)    # config4, all datasets summary


exon.list = list()
intron.list = list()
df.fisher = NULL
df.dynTFs = NULL
df.pval = data.frame()
par(mfrow =c(1,1))

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

# # load ground-truth
if(T){
  # gt <--> sample
  if(sample_source=="human"){
    gt_id = 6
  }
  if(sample_source=="mouse"){
    gt_id = 7
  }
  
  # ground-truth network
  if(T){
    gt_names = c("dorothea","RegM","TRRUST_M","Doro_M","motif_sel","motif","motif_m","own_gt")
    gt_name = gt_names[gt_id]
    #ground-truth -- dorothea
    if(gt_id == 1){
      load(file ="./dorothea/dorothea_database.RData")
      tmp = unique(dorothea_hs_high_confidence)
      tmp = as.matrix(tmp)
      reference = cbind(as.character(tmp[,1]),as.character(tmp[,3]))
      gt = list_to_mat(reference)
    }
    if(gt_id == 2){
      load(file = "./RegNetwork/RegNetwork_mat_mouse.RData")
      gt = mat
    }
    if(gt_id == 3){
      load(file = "./TRRUST/TRRUST_mouse.RData")
      gt = mat
    }
    if(gt_id == 4){
      load(file ="./dorothea/dorothea_database.RData")
      tmp = unique(dorothea_mm_high_confidence)
      tmp = as.matrix(tmp)
      reference = cbind(as.character(tmp[,1]),as.character(tmp[,3]))
      gt = list_to_mat(reference)
    }
    # motif selected network
    if(gt_id==5){
      load(file = "./hESC_perturb/motif_selected_network.RData")
      gt = motif_sel_network
    }
    # motif network
    if(gt_id==6){
      load(file = "./Rcistarget/motif_link_and_mat.RData")
      gt = motif_mat
    }
    # motif network -- mouse
    if(gt_id==7){
      load(file = "./Rcistarget/mESC/motif_link_and_mat.RData")
      gt = motif_mat
    }
    # own gt: ChIP X RNA-seq
    if(gt_id==8){
      load(file = "./construct_GRN/0731_chipxRNAseq_net.RData")
      gt = mat_list[[own_id]]
    }
    
    
    # remove autoregulation
    TFs_col_id = match(rownames(gt),colnames(gt))
    sub_id = which(!is.na(TFs_col_id))
    for(i in sub_id){
      gt[i,TFs_col_id[i]] = 0
    }
    
  }
  
  # change size
  if(T){
    dim(exon.mat)
    dim(inex.mat)
    overlap_targets = intersect(colnames(exon.mat),colnames(inex.mat))
    exon.mat = exon.mat[,overlap_targets]
    inex.mat = inex.mat[,overlap_targets]
    dim(inex.mat)
    
    dim(gt)
    TFs = intersect(rownames(gt),rownames(inex.mat))
    targets = intersect(colnames(gt),colnames(inex.mat))
    gt = gt[TFs,targets]
    dim(gt)
    sum(as.vector(gt))
    exon.mat = exon.mat[TFs,targets]
    inex.mat = inex.mat[TFs,targets]
    
  }
}

# load mRNA count number data
if(T){
  library(Matrix)
  if(sample_id==1){
    dir = './10X_datasets/unspliced_spliced_matrix/Velocyto/Forebrain'
  }
  if(sample_id==2){
    dir = './10X_datasets/unspliced_spliced_matrix/Velocyto/BoneMarrow'
  }
  if(sample_id==3){
    dir = './10X_datasets/unspliced_spliced_matrix/Velocyto/epithelium'
  }
  if(sample_id==4){
    dir = './10X_datasets/unspliced_spliced_matrix/scVelo/DentateGyrus'
  }
  if(sample_id==5){
    dir = './10X_datasets/unspliced_spliced_matrix/scVelo/endocrinogenesis'
  }
  if(sample_id==13){
    dir = './RNAseq/hESC/SRR6328624_2018'
  }
  if(sample_id==14){
    dir = './RNAseq/hESC/SRR9117953_2019'
  }
  if(sample_id==15){
    dir = './RNAseq/GSE152125'
  }
  # load raw data
  load(file = paste(dir,'/unspliced_spliced_matrix.RData',sep = ''))
  u.raw = u.mat
  s.raw = s.mat
  # mean count
  u.m = rowSums(u.raw)/ncol(u.raw)
  s.m = rowSums(s.raw)/ncol(s.raw)
  u.sd = apply(u.raw,1,sd)
  s.sd = apply(s.raw,1,sd)
  # mapping gene names
  if(is.element(sample_id,c(13,14))){
    library(AnnotationDbi)
    library(org.Hs.eg.db) 
    geneset = mapIds(org.Hs.eg.db,names(u.m),keytype = "ENSEMBL",column = "SYMBOL")
    length(geneset)
    length(na.omit(geneset))
    id = which(!is.na(geneset))
    # change names
    u.m = u.m[id]
    s.m = s.m[id]
    u.sd = u.sd[id]
    s.sd = s.sd[id]
    names(u.m) = names(s.m) = names(u.sd) = names(s.sd) = geneset[id]
  }
  if(sample_id ==15){
    names(u.m) = names(u.sd) = sub("-intron","",names(u.m))
  }
}
# evaluation with each TF
if(F){
  
  # Evaluation
  # df.exon = Evaluation_TF(exon.mat,gt)
  # df.inex = Evaluation_TF(inex.mat,gt)
  
  # top_num = 50
  # tmp = Eval_top_for_TFs(exon.mat,gt,top_num)
  # df.exon$top_precision = tmp$top_precision
  # tmp = Eval_top_for_TFs(inex.mat,gt,top_num)
  # df.inex$top_precision = tmp$top_precision
  
  # save(df.exon,df.inex,file = "./Evaluation_with_hub_TFs/evaluation_results/hFB_motif_eval_TFs.RData")
  # save(df.exon,df.inex,file = "./Evaluation_with_hub_TFs/evaluation_results/mFB_motif_eval_TFs.RData")
  filename = paste("./Evaluation_with_hub_TFs/evaluation_results/",sample_name,"_motif_eval_TFs.RData",sep = "")
  # save(df.exon,df.inex,file = filename)
  
  load(file = filename)
  
  # id.zero = which(df.exon$n_true_targets==0)
  # df.exon = df.exon[-id.zero,]
  # df.inex = df.inex[-id.zero,]
  
  # investigate the effect of mean count
  if(F){
    # map to subset
    match_id = match(rownames(gt),names(u.m))
    u.m.sub = u.m[match_id]
    match_id = match(rownames(gt),names(s.m))
    s.m.sub = s.m[match_id]
    
    # meancount - better proportion regression
    df.inex.all = df.inex
    df.inex.all$m = s.m.sub[match(rownames(df.inex.all),names(s.m.sub))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    
    # calculate proportion
    cal_TP = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.tp = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$m>values[i]&df.inex.all$m<values[i+1])
        tp = mean(df.inex.all$top_precision[id])
        diff = mean(df.inex.all$diff[id])
        random_precision = mean(df.inex.all$random_precision[id])
        tmp = c(values[i],values[i+1],length(id),tp,random_precision,diff)
        df.tp = rbind(df.tp,tmp)
      }
      colnames(df.tp) = c("value_min","value_max","n","TP","random_precision","diff")
      df.tp = as.data.frame(df.tp)
      return(df.tp)
    }
    
    # divide into groups to calculate proportion
    n_bins = 10
    values = 2**seq(log2(min(df.inex.all$m)),log2(max(df.inex.all$m+1)),length.out = n_bins)
    
    df.tp.in = cal_TP(df.inex.all,values)
    # exon results
    df.exon.all = df.exon
    df.exon.all$m = s.m.sub[match(rownames(df.exon.all),names(s.m.sub))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    df.tp.ex = cal_TP(df.exon.all,values)
    
    plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$TP,type = "b",col = "red",log = "x",
         xlim = range(values),ylim = c(0,max(c(df.tp.in$TP,df.tp.ex$TP)+0.05,na.rm = T)),xlab = "mean count", ylab = "Top precision",
         main = paste("Effect of mean count on TF: ",sample_name,sep = ""))
    
    
    
    lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$TP,type = "b",col = "blue")
    lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$random_precision,col = "black",lty = 2)
    
    
  }
  # investigate the effect of sd
  if(F){
    # map to subset
    match_id = match(rownames(gt),names(s.sd))
    s.sd.sub = s.sd[match_id]
    
    # meancount - better proportion regression
    df.inex.all = df.inex
    df.inex.all$sd = s.sd.sub[match(rownames(df.inex.all),names(s.sd.sub))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    
    # calculate proportion
    cal_TP = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.tp = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$sd>values[i]&df.inex.all$sd<values[i+1])
        tp = mean(df.inex.all$top_precision[id])
        diff = mean(df.inex.all$diff[id])
        random_precision = mean(df.inex.all$random_precision[id])
        tmp = c(values[i],values[i+1],length(id),tp,random_precision,diff)
        df.tp = rbind(df.tp,tmp)
      }
      colnames(df.tp) = c("value_min","value_max","n","TP","random_precision","diff")
      df.tp = as.data.frame(df.tp)
      return(df.tp)
    }
    # divide into groups to calculate proportion
    values = 2**seq(log2(min(df.inex.all$sd)),log2(max(df.inex.all$sd+1)),length.out = n_bins)
    
    df.tp.in = cal_TP(df.inex.all,values)
    plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$TP,type = "b",col = "red",log = "x",
         xlim = range(values),ylim = c(0,max(c(df.tp.in$TP,df.tp.ex$TP)+0.05,na.rm = T)),xlab = "sd", ylab = "Top precision",
         main = paste("Effect of sd on TF: ",sample_name,sep = ""))
    
    # exon results
    df.exon.all = df.exon
    df.exon.all$sd = s.sd.sub[match(rownames(df.exon.all),names(s.sd.sub))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision

    df.tp.ex = cal_TP(df.exon.all,values)
    lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$TP,type = "b",col = "blue")
    lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$random_precision,col = "black",lty = 2)
  }
  # investigate the effect of variability
  if(F){
    # load data
    load(file = "./ATACseq_var/processed_vector.RData")
    var.m = var_m_TF
    overlap_genes = intersect(toupper(rownames(df.exon)),names(var.m))
    # meancount - better proportion regression
    df.inex.all = df.inex[overlap_genes,]
    df.inex.all$var = var.m[overlap_genes]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    
    # calculate proportion
    cal_TP = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.tp = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$var>values[i]&df.inex.all$var<values[i+1])
        tp = mean(df.inex.all$top_precision[id])
        diff = mean(df.inex.all$diff[id])
        random_precision = mean(df.inex.all$random_precision[id])
        tmp = c(values[i],values[i+1],length(id),tp,random_precision,diff)
        df.tp = rbind(df.tp,tmp)
      }
      colnames(df.tp) = c("value_min","value_max","n","TP","random_precision","diff")
      df.tp = as.data.frame(df.tp)
      return(df.tp)
    }
    # divide into groups to calculate proportion
    n_bins = 10
    values = 2**seq(log2(min(df.inex.all$var)),log2(max(df.inex.all$var)),length.out = n_bins)
    
    df.tp.in = cal_TP(df.inex.all,values)
    
    # exon results
    df.exon.all = df.exon[overlap_genes,]
    df.exon.all$var = var.m[overlap_genes]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    
    df.tp.ex = cal_TP(df.exon.all,values)
    
    plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$TP,type = "b",col = "red",log = "x",
         xlim = range(values),ylim = c(0,max(c(df.tp.in$TP,df.tp.ex$TP)+0.05,na.rm = T)),xlab = "variability", 
         ylab = "Top precision",main = paste("Effect of variability on TF: ",sample_name,sep = ""))
    lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$TP,type = "b",col = "blue")
    lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$random_precision,col = "black",lty = 2)
  }
  
  # investigate dynamic TFs
  if(T){
    df.diff = df.inex
    df.diff$diff = df.inex$AUPR - df.exon$AUPR
    # focus on at least one right TF
    id = which(df.inex$AUPR<df.inex$random_precision&df.exon$AUPR<df.inex$random_precision)
    df.diff = df.diff[-id,]
    # # select hub TFs
    # df.diff = df.diff[intersect(hub_TFs,rownames(df.diff)),]
    rownames(df.diff) = toupper(rownames(df.diff))
    
    # load data
    if(F){
      load(file = "./ATACseq_var/processed_matrix.RData")
      var.m = apply(var_mat,1,mean)
    }
    load(file = "./ATACseq_var/processed_vector.RData")
    var.m = var_m_TF
    # load variability for human forebrain
    if(F){
      load(file = "./ATACseq_var/forebrain/processed_var_vec.RData")
      var.m = var_vec
      # hist(var.m)
    }
    overlap_genes = intersect(toupper(rownames(df.diff)),names(var.m))
    
    tmp = rep(NA,length(var.m))
    tmp[match(overlap_genes,names(var.m))] = df.diff[overlap_genes,"diff"]
    df.dynTFs = cbind(df.dynTFs,tmp)
    
    df.diff = df.diff[overlap_genes,]
    df.diff$variability = var.m[rownames(df.diff)]
    
    intron_TFs = rownames(df.diff)[which(df.diff$diff>0)]
    exon_TFs = rownames(df.diff)[which(df.diff$diff<0)]
    
    
    
    # par(mfrow = c(1,2))
    # plot(df.diff$diff,df.diff$variability,main = paste("Intron vs Exon",sample_name,sep = " "),xlab = "AUPR: Intron - exon",
    #      ylab = "variability")
    # abline(v = 0,lty = 2)
    # abline(h = 1.2,lty = 2)
    # 
    # boxplot(list(a = var.m[intron_TFs],b = var.m[exon_TFs]),names = c("intron_TFs","exon_TFs"),ylab = "variability",
    #         main = paste("Intron vs Exon",sample_name,sep = " "))
    pvalue = wilcox.test(var.m[intron_TFs],var.m[exon_TFs])$p.value
    print(paste("pvalue: ",pvalue,sep = ""))
    
    print(paste(sample_name,", Intron better variable TFs: ",sep = ""))
    print(rownames(df.diff)[which(df.diff$diff>0&df.diff$variability>1.2)])
    print(paste(sample_name,", Exon better variable TFs: ",sep = ""))
    print(rownames(df.diff)[which(df.diff$diff<0&df.diff$variability>1.2)])
    
    n_in_nonvar = length(which(df.diff$diff>0&df.diff$variability<=1.2))
    n_ex_nonvar = length(which(df.diff$diff<0&df.diff$variability<=1.2))
    n_in_var = length(which(df.diff$diff>0&df.diff$variability>1.2))
    n_ex_var = length(which(df.diff$diff<0&df.diff$variability>1.2))
    a = fisher.test(rbind(c(n_ex_nonvar,n_in_nonvar),c(n_ex_var,n_in_var)))
    odds_ratio = a$estimate
    pvalue = a$p.value
    df.fisher = rbind(df.fisher,c(n_ex_nonvar,n_in_nonvar,n_ex_var,n_in_var,odds_ratio,pvalue))
    
    
    
    # plot proportion - variability 
    if(T){
      tmp = rep(NA,length(var.m))
      tmp[match(overlap_genes,names(var.m))] = df.diff[overlap_genes,"diff"]
      df.tmp = cbind(tmp,var.m)
      rownames(df.tmp) = names(var.m)
      colnames(df.tmp) = c("diff","var")
      df.tmp = as.data.frame(df.tmp)
      df.tmp = df.tmp[order(df.tmp$var,decreasing = T),]
      # calculate number of sample that have overlap
      n = apply(df.tmp,1,function(x){return(1 - sum(is.na(x)))})
      df.tmp$n =n
      # calculate intron better number and exon better number
      n_in = sapply(df.tmp[,1],function(x){return(length(which(x>0)))})
      n_ex = sapply(df.tmp[,1],function(x){return(length(which(x<0)))})
      df.tmp$n_in =n_in
      df.tmp$n_ex =n_ex
      
      # calculate variability - proportion results
      df.sub = df.tmp[which(df.tmp$n>0),]
      n_in_sum = n_ex_sum = 0
      props = NULL
      for(i in 1:nrow(df.sub)){
        n_in_sum = n_in_sum + df.sub$n_in[i]
        n_ex_sum = n_ex_sum + df.sub$n_ex[i]
        props = c(props,n_in_sum/(n_in_sum+n_ex_sum))
      }
      plot(df.sub$var,props,type = "b",main = paste("Dynamic TFs",sample_names[sample_id]),xlab = "variability",
           ylab = "proportion: AUPR intron > exon",ylim= range(c(0.5,props),na.rm = T))
      abline(h = 0.5,col = "grey")
    }
    
    
    # summary 4 datasets
    if(F){
      colnames(df.fisher) = c("n_ex_nonvar","n_in_nonvar","n_ex_var","n_in_var","odds_ratio","pvalue")
      rownames(df.fisher) = sample_names[sample_ids]
      df.fisher= as.data.frame(df.fisher)
      
      n_ex_nonvar = sum(df.fisher$n_ex_nonvar)
      n_in_nonvar = sum(df.fisher$n_in_nonvar)
      n_ex_var = sum(df.fisher$n_ex_var)
      n_in_var = sum(df.fisher$n_in_var)
      a = fisher.test(rbind(c(n_ex_nonvar,n_in_nonvar),c(n_ex_var,n_in_var)))
      odds_ratio = a$estimate
      pvalue = a$p.value
      
      
      # evaluate dynamic TFs
      rownames(df.dynTFs) = names(var.m)
      colnames(df.dynTFs) = sample_names[sample_ids]
      df.dynTFs = as.data.frame(df.dynTFs)
      df.dynTFs$var = var.m
      df.dynTFs = df.dynTFs[order(df.dynTFs$var,decreasing = T),]
      # calculate number of sample that have overlap
      n = apply(df.dynTFs,1,function(x){return(4 - sum(is.na(x)))})
      df.dynTFs$n =n
      id = which(df.dynTFs$n>=1)
      write.csv(df.dynTFs[id,1:4],file = "TF_variability_4_datasets.csv",quote = F)
      View(df.dynTFs[id,])
      # calculate intron better number and exon better number
      n_in = apply(df.dynTFs[,1:4],1,function(x){return(length(which(x>0)))})
      n_ex = apply(df.dynTFs[,1:4],1,function(x){return(length(which(x<0)))})
      df.dynTFs$n_in =n_in
      df.dynTFs$n_ex =n_ex
      
      # calculate variability - proportion results
      df.sub = df.dynTFs[which(df.dynTFs$n>0),]
      n_in_sum = n_ex_sum = 0
      props = NULL
      for(i in 1:nrow(df.sub)){
        n_in_sum = n_in_sum + df.sub$n_in[i]
        n_ex_sum = n_ex_sum + df.sub$n_ex[i]
        props = c(props,n_in_sum/(n_in_sum+n_ex_sum))
      }
      plot(df.sub$var,props,type = "b",main = "Dynamic TFs show intron better than exon",xlab = "variability",
           ylab = "proportion: AUPR intron > exon")
      abline(h = 0.5,col = "grey")
      
      
      # calculate barplot
      if(F){
        thr = 1.2
        id.var = which(df.dynTFs$var>thr)
        id.nonvar = which(df.dynTFs$var<=thr)
        
        id = id.var
        n_in = sum(df.dynTFs$n_in[id])
        n_ex = sum(df.dynTFs$n_ex[id])
        n = n_in + n_ex
        prop = n_in/n
        df_var = data.frame(n = n, prop = prop)
        
        id = id.nonvar
        n_in = sum(df.dynTFs$n_in[id])
        n_ex = sum(df.dynTFs$n_ex[id])
        n = n_in + n_ex
        prop = n_in/n
        df_nonvar = data.frame(n = n, prop = prop)
        
        barplot(c(df_var$prop,df_nonvar$prop),names.arg = c("variable TFs","non-variable TFs"),col = "white",
                ylab = "proportion: intron better than exon",main = "Compare intron with exon for dynamic TFs",
                ylim = c(0,0.8))
        abline(h = 0.5,col ="grey",lty = 2)
        
        # boxplot for AUPR difference
        diff.var = na.omit(as.vector(as.matrix(df.dynTFs[id.var,1:4])))
        diff.nonvar = na.omit(as.vector(as.matrix(df.dynTFs[id.nonvar,1:4])))
        boxplot(list(a = diff.var,b = diff.nonvar),names = c("variable TFs","non-variable TFs"),
                ylab = "AUPR difference: intron - exon",main = "Compare intron with exon for dynamic TFs")
        abline(h = 0,col ="grey",lty = 2)
        # t.test(diff.var)
        wilcox.test(diff.var,diff.nonvar,alternative = "greater")
      }
    }
  }
  
}

# evaluation with each target
if(T){
  # Evaluation
  # df.exon = Evaluation_target(exon.mat,gt)
  # df.inex = Evaluation_target(inex.mat,gt)
  
  # top_num = 10
  # tmp = Eval_top_for_targets(exon.mat,gt,top_num)
  # df.exon$top_precision = tmp$top_precision
  # tmp = Eval_top_for_targets(inex.mat,gt,top_num)
  # df.inex$top_precision = tmp$top_precision
  
  # save(df.exon,df.inex,file = "./Evaluation_with_hub_TFs/evaluation_results/mFB_motif_eval_targets.RData")
  filename = paste("./Evaluation_with_hub_TFs/evaluation_results/",sample_name,"_motif_eval_targets.RData",sep = "")
  # save(df.exon,df.inex,file = filename)
  
  load(file = filename)
  
  # explore the effect of mean count
  if(T){
    # map to subset
    match_id = match(colnames(gt),names(u.m))
    u.m.sub = u.m[match_id]
    match_id = match(colnames(gt),names(s.m))
    s.m.sub = s.m[match_id]
    
    # meancount - better proportion regression
    df.inex.all = df.inex
    df.inex.all$m = u.m.sub[match(rownames(df.inex.all),names(u.m.sub))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    df.inex.all = df.inex.all[-which(df.inex.all$n_true_TFs==0),]
    
    # plot(df.inex.all$m,df.inex.all$diff,log = "x")
    # plot(df.inex.all$m,df.inex.all$top_precision,log = "x")
    
    
    # calculate proportion
    cal_prop = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.prop = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$m>values[i]&df.inex.all$m<values[i+1])
        better.id = intersect(which(df.inex.all$diff>0),id)
        worse.id = intersect(which(df.inex.all$diff<0),id)
        prop = length(better.id)/(length(better.id)+length(worse.id))
        tmp = c(values[i],values[i+1],length(better.id),length(worse.id),prop)
        df.prop = rbind(df.prop,tmp)
      }
      colnames(df.prop) = c("value_min","value_max","n_better","n_worse","prop")
      df.prop = as.data.frame(df.prop)
      return(df.prop)
    }
    cal_TP = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.tp = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$m>values[i]&df.inex.all$m<values[i+1])
        tp = mean(df.inex.all$top_precision[id])
        random_precision = mean(df.inex.all$random_precision[id])
        diff = mean(df.inex.all$diff[id],na.rm = T)
        tmp = c(values[i],values[i+1],length(id),tp,random_precision,diff)
        df.tp = rbind(df.tp,tmp)
      }
      colnames(df.tp) = c("value_min","value_max","n","TP","random_precision","diff")
      df.tp = as.data.frame(df.tp)
      return(df.tp)
    }
    
    
    # exon results
    df.exon.all = df.exon
    df.exon.all$m = s.m.sub[match(rownames(df.exon.all),names(s.m.sub))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    df.exon.all = df.exon.all[-which(df.exon.all$n_true_TFs==0),]
    
    # plot(df.exon.all$m,df.exon.all$diff>0,log = "x")
    
    # divide into groups to calculate proportion
    if(F){
      n_bins = 10
      values = c(0.1,0.2,0.5,1,2,5,10,20)
      values = 2**seq(log2(min(df.inex.all$m)),log2(max(df.inex.all$m+1)),length.out = n_bins)
      df.prop.in = cal_prop(df.inex.all,values)
      plot((df.prop.in$value_min+df.prop.in$value_max)/2,df.prop.in$prop,type = "b",col = "red",log = "x",
           xlim = c(min(0.1,df.inex.all$m),300),ylim = c(0,1),xlab = "mean count", ylab = "Better than random proportion",
           main = "Effect of mean count for inference")
      
      values = c(0.1,0.2,0.5,1,2,5,10,20,50,100)
      values = 2**seq(log2(min(df.exon.all$m)),log2(max(df.exon.all$m+1)),length.out = n_bins)
      
      df.prop.ex = cal_prop(df.exon.all,values)
      lines((df.prop.ex$value_min+df.prop.ex$value_max)/2,df.prop.ex$prop,type = "b",col = "blue")
    }
    
    # divide into groups to calculate Top precision
    if(T){
      n_bins = 10
      
      values = 2**seq(log2(min(df.inex.all$m)),log2(max(df.inex.all$m+1)),length.out = n_bins)
      
      df.tp.in = cal_TP(df.inex.all,values)
      values = 2**seq(log2(min(df.exon.all$m)),log2(max(df.exon.all$m+1)),length.out = n_bins)
      
      df.tp.ex = cal_TP(df.exon.all,values)
      plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$TP,type = "b",col = "red",log = "x",
           xlim = c(min(0.1,df.inex.all$m),300),ylim = c(0,max(c(df.tp.in$TP,df.tp.ex$TP)+0.05,na.rm = T)),xlab = "mean count", ylab = "Top precision",
           main = paste("Effect of mean count: ",sample_name,sep = ""))
      
      
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$TP,type = "b",col = "blue")
      lines((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$random_precision,type = "l",col = "red",lty = 2)
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$random_precision,type = "l",col = "blue",lty = 2)
      
      
    }
    
    # divide into groups to calculate AUPR difference
    if(F){
      n_bins = 10
      
      values = 2**seq(log2(min(df.inex.all$m)),log2(max(df.inex.all$m+1)),length.out = n_bins)
      
      df.tp.in = cal_TP(df.inex.all,values)
      values = 2**seq(log2(min(df.exon.all$m)),log2(max(df.exon.all$m+1)),length.out = n_bins)
      
      df.tp.ex = cal_TP(df.exon.all,values)
      plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$diff,type = "b",col = "red",log = "x",
           xlim = c(min(0.1,df.inex.all$m),300),ylim = range(c(df.tp.in$diff,df.tp.ex$diff,0),na.rm = T),
           xlab = "mean count", ylab = "AUPR difference",
           main = paste("Effect of mean count: ",sample_name,sep = ""))
      
      
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$diff,type = "b",col = "blue")
      abline(h = 0,lty = 2,col = "grey")
      
      
    }
    
    
  }
  
  # explore the effect of mean count: AUROC
  if(T){
    # map to subset
    match_id = match(colnames(gt),names(u.m))
    u.m.sub = u.m[match_id]
    match_id = match(colnames(gt),names(s.m))
    s.m.sub = s.m[match_id]
    
    # meancount - better proportion regression
    df.inex.all = df.inex
    df.inex.all$m = u.m.sub[match(rownames(df.inex.all),names(u.m.sub))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    df.inex.all = df.inex.all[-which(df.inex.all$n_true_TFs==0),]
    
    # plot(df.inex.all$m,df.inex.all$diff,log = "x")
    # plot(df.inex.all$m,df.inex.all$top_precision,log = "x")
    
    
    # calculate proportion
    cal_prop = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.prop = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$m>values[i]&df.inex.all$m<values[i+1])
        better.id = intersect(which(df.inex.all$diff>0),id)
        worse.id = intersect(which(df.inex.all$diff<0),id)
        prop = length(better.id)/(length(better.id)+length(worse.id))
        tmp = c(values[i],values[i+1],length(better.id),length(worse.id),prop)
        df.prop = rbind(df.prop,tmp)
      }
      colnames(df.prop) = c("value_min","value_max","n_better","n_worse","prop")
      df.prop = as.data.frame(df.prop)
      return(df.prop)
    }
    cal_TP = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.tp = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$m>values[i]&df.inex.all$m<values[i+1])
        tp = mean(df.inex.all$top_precision[id])
        random_precision = mean(df.inex.all$random_precision[id])
        diff = mean(df.inex.all$diff[id],na.rm = T)
        tmp = c(values[i],values[i+1],length(id),tp,random_precision,diff)
        df.tp = rbind(df.tp,tmp)
      }
      colnames(df.tp) = c("value_min","value_max","n","TP","random_precision","diff")
      df.tp = as.data.frame(df.tp)
      return(df.tp)
    }
    cal_auroc = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.auroc = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$m>values[i]&df.inex.all$m<values[i+1])
        AUROC = mean(df.inex.all$AUROC[id])
        random_precision = mean(df.inex.all$random_precision[id])
        diff = mean(df.inex.all$diff[id],na.rm = T)
        tmp = c(values[i],values[i+1],length(id),AUROC,random_precision,diff)
        df.auroc = rbind(df.auroc,tmp)
      }
      colnames(df.auroc) = c("value_min","value_max","n","AUROC","random_precision","diff")
      df.auroc = as.data.frame(df.auroc)
      return(df.auroc)
    }
    
    
    # exon results
    df.exon.all = df.exon
    df.exon.all$m = s.m.sub[match(rownames(df.exon.all),names(s.m.sub))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    df.exon.all = df.exon.all[-which(df.exon.all$n_true_TFs==0),]
    
    # plot(df.exon.all$m,df.exon.all$diff>0,log = "x")
    
    # divide into groups to calculate AUROC
    if(T){
      n_bins = 10
      
      values = 2**seq(log2(min(df.inex.all$m)),log2(max(df.inex.all$m+1)),length.out = n_bins)
      
      df.tp.in = cal_auroc(df.inex.all,values)
      values = 2**seq(log2(min(df.exon.all$m)),log2(max(df.exon.all$m+1)),length.out = n_bins)
      
      df.tp.ex = cal_auroc(df.exon.all,values)
      plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$AUROC,type = "b",col = "red",log = "x",
           xlim = c(min(0.1,df.inex.all$m),300),ylim = c(0,max(c(df.tp.in$AUROC,df.tp.ex$AUROC)+0.05,na.rm = T)),xlab = "mean count", ylab = "AUROC",
           main = paste("Effect of mean count: ",sample_name,sep = ""))
      
      
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$AUROC,type = "b",col = "blue")
      
      
    }
    
    
  }
  
  # load mRNA half-life data
  if(T){
    # K562 data
    if(T){
      a = openxlsx::read.xlsx("./kinetic_parameters/mRNA_hl/2018 time-lapse-seq.xlsx",sheet =1)
      hl = a$mean_half_life
      names(hl) = a$transcript
    }
    # match gene names
    rownames(df.inex) = toupper(rownames(df.inex))
    rownames(df.exon) = toupper(rownames(df.exon))
    rownames(gt) = toupper(rownames(gt))
    colnames(gt) = toupper(colnames(gt))
    
    dim(gt)
    match_id = match(rownames(df.exon),names(hl))
    # match_id = match(colnames(gt),names(hl))
    
    hl.matched = hl[match_id]
    length(hl.matched)
    length(na.omit(hl.matched))
    # genes taht have mRNA half-life data
    id.hldata = which(!is.na(hl.matched))
    
  }
  
  # effect of mRNA half-life -- detailed view
  if(T){
    # mRNA half-life - better proportion regression
    df.inex.all = df.inex
    df.inex.all$hl = hl.matched[match(rownames(df.inex.all),names(hl.matched))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision

    nrow(df.inex.all)
    nrow(na.omit(df.inex.all))
    
    # plot(df.inex.all$hl,df.inex.all$diff)
    # calculate proportion
    cal_prop = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.prop = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$hl>values[i]&df.inex.all$hl<values[i+1])
        better.id = intersect(which(df.inex.all$diff>0),id)
        worse.id = intersect(which(df.inex.all$diff<0),id)
        prop = length(better.id)/(length(better.id)+length(worse.id))
        tmp = c(values[i],values[i+1],length(better.id),length(worse.id),prop)
        df.prop = rbind(df.prop,tmp)
      }
      colnames(df.prop) = c("value_min","value_max","n_better","n_worse","prop")
      df.prop = as.data.frame(df.prop)
      return(df.prop)
    }
    cal_prop_inex = function(df.inex.all,df.exon.all,values){
      n_groups = length(values)-1
      df.prop = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$hl>values[i]&df.inex.all$hl<values[i+1])
        better.id = intersect(which(df.inex.all$AUPR - df.exon.all$AUPR>0),id)
        worse.id = intersect(which(df.inex.all$AUPR - df.exon.all$AUPR<0),id)
        prop = length(better.id)/(length(better.id)+length(worse.id))
        tmp = c(values[i],values[i+1],length(better.id),length(worse.id),prop)
        df.prop = rbind(df.prop,tmp)
      }
      colnames(df.prop) = c("value_min","value_max","n_better","n_worse","prop")
      df.prop = as.data.frame(df.prop)
      return(df.prop)
    }
    cal_TP = function(df,values){
      df.inex.all = df
      n_groups = length(values)-1
      df.tp = NULL
      for(i in 1:n_groups){
        id = which(df.inex.all$hl>values[i]&df.inex.all$hl<values[i+1])
        tp = mean(df.inex.all$top_precision[id])
        diff = mean(df.inex.all$diff[id],na.rm = T)
        random_precision = mean(df.inex.all$random_precision[id])
        tmp = c(values[i],values[i+1],length(id),tp,random_precision,diff)
        df.tp = rbind(df.tp,tmp)
      }
      colnames(df.tp) = c("value_min","value_max","n","TP","random_precision","diff")
      df.tp = as.data.frame(df.tp)
      return(df.tp)
    }
    
    # exon results
    df.exon.all = df.exon
    df.exon.all$hl = hl.matched[match(rownames(df.exon.all),names(hl.matched))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    
    
    
    # divide into groups to calculate proportion
    if(F){
      nbins = 10
      values = 2**seq(log2(min(na.omit(df.inex.all$hl))),log2(max(na.omit(df.inex.all$hl))+0.1),length.out = n_bins)
      
      df.prop.in = cal_prop(df.inex.all,values)
      plot((df.prop.in$value_min+df.prop.in$value_max)/2,df.prop.in$prop,type = "b",col = "red",log = "x",
           xlim = c(0.2,20),ylim = c(0,1),xlab = "mRNA half-life (h)", ylab = "Better than random proportion",
           main = "Effect of mRNA half-life for inference")
      
      df.prop.ex = cal_prop(df.exon.all,values)
      lines((df.prop.ex$value_min+df.prop.ex$value_max)/2,df.prop.ex$prop,type = "b",col = "blue")
    }
    
    # divide into groups to calculate Top precision
    if(T){
      nbins = 10
      values = 2**seq(log2(min(na.omit(df.inex.all$hl))),log2(max(na.omit(df.inex.all$hl))+0.1),length.out = n_bins)
      
      df.tp.in = cal_TP(df.inex.all,values)
      df.tp.ex = cal_TP(df.exon.all,values)
      
      plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$TP,type = "b",col = "red",log = "x",
           xlim = c(0.2,20),ylim = c(0,max(c(df.tp.in$TP,df.tp.ex$TP)+0.05,na.rm = T)),xlab = "mRNA half-life (h)", ylab = "Top precision",
           main = paste("Effect of mRNA half-life: ",sample_name,sep = ""))
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$TP,type = "b",col = "blue")
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$random_precision,type = "l",col = "black",lty = 2)
      
    }
    
    # divide into groups to calculate AUPR difference
    if(F){
      n_bins = 10
      
      values = 2**seq(log2(min(na.omit(df.inex.all$hl))),log2(max(na.omit(df.inex.all$hl))+0.1),length.out = n_bins)
      
      df.tp.in = cal_TP(df.inex.all,values)
      df.tp.ex = cal_TP(df.exon.all,values)
      
      plot((df.tp.in$value_min+df.tp.in$value_max)/2,df.tp.in$diff,type = "b",col = "red",log = "x",
           xlim = c(min(0.1,df.inex.all$m),300),ylim = range(c(df.tp.in$diff,df.tp.ex$diff,0),na.rm = T),
           xlab = "mRNA half-life (h)", ylab = "AUPR difference",
           main = paste("Effect of mean count: ",sample_name,sep = ""))
      
      
      lines((df.tp.ex$value_min+df.tp.ex$value_max)/2,df.tp.ex$diff,type = "b",col = "blue")
      abline(h = 0,lty = 2,col = "grey")
      
      
    }
  }
  
  # correlation of meancount and half-life
  if(F){
    names(u.m) = toupper(names(u.m))
    names(s.m) = toupper(names(s.m))
    # map to subset
    match_id = match(colnames(gt),names(u.m))
    u.m.sub = u.m[match_id]
    match_id = match(colnames(gt),names(s.m))
    s.m.sub = s.m[match_id]
    
    hl.tmp = hl.matched[id.hldata]
    s.m.tmp = s.m.sub[id.hldata]
    u.m.tmp = u.m.sub[id.hldata]
    
    id = which(!is.na(s.m.tmp)&!is.na(u.m.tmp))
    hl.tmp = hl.tmp[id]
    s.m.tmp = s.m.tmp[id]
    u.m.tmp = u.m.tmp[id]
    
    par(mfrow = c(2,2))
    plot(hl.tmp,s.m.tmp,col ="blue",main = paste("cor =",round(cor(hl.tmp,log2(s.m.tmp)),2),sample_name,sep = " "),log = "y",
         xlab = "mRNA half-life",ylab = "mean count")
    plot(hl.tmp,u.m.tmp,col ="red",main = paste("cor =",round(cor(hl.tmp,log2(u.m.tmp)),2),sample_name,sep = " "),log = "y",
         xlab = "mRNA half-life",ylab = "mean count")
    # plot(s.m.tmp,u.m.tmp,col ="red",main = paste("cor =",round(cor(s.m.tmp,u.m.tmp),2),sample_name,sep = " "),log = "xy")
    # cor(log2(s.m.tmp),log2(u.m.tmp))
    cor(hl.tmp,s.m.tmp)
    cor(hl.tmp,s.m.tmp,method = "spearman")
    
    
    df.targets = data.frame(exon.m = s.m.tmp,inex.m = u.m.tmp, hl = hl.tmp)
    
    values = seq(0,24,2)
    df = NULL
    a = b = list()
    for(i in 1:(length(values)-1)){
      id = which(df.targets$hl>values[i]&df.targets$hl<=values[i+1])
      a[[i]] = df.targets$exon.m[id]
      b[[i]] = df.targets$inex.m[id]
      df = rbind(df,c(length(id),values[i],values[i+1],mean(df.targets$exon.m[id]),mean(df.targets$inex.m[id])))
    }
    colnames(df) = c("n_targets","min_value","max_value","exon.m","inex.m")
    df = as.data.frame(df)
    boxplot(a,log = "y")
    boxplot(b,log = "y")
    
    # plot((df$min_value+df$max_value)/2,df$exon.m,main = paste("exon.m - mRNA half-life: ",sample_name,sep = ""),
    #      xlab = "mRNA half-life",ylab = "mean count",type = "b",xlim = c(0,24),log = "y")
    # plot((df$min_value+df$max_value)/2,df$inex.m,main = paste("inex.m - mRNA half-life: ",sample_name,sep = ""),
    #      xlab = "mRNA half-life",ylab = "mean count",type = "b",xlim = c(0,24),log = "y")
    # 
    # wilcox.test(c(b[[1]],b[[2]]),unlist(b))
  }
  # investigate the effect of mean count and mRNA half-life
  if(F){
    par(mfrow = c(1,4))
    # match gene names
    rownames(df.inex) = toupper(rownames(df.inex))
    rownames(df.exon) = toupper(rownames(df.exon))
    rownames(gt) = toupper(rownames(gt))
    colnames(gt) = toupper(colnames(gt))
    names(s.m) = toupper(names(s.m))
    names(u.m) = toupper(names(u.m))
    
    # mRNA half-life - better proportion regression
    df.exon.all = df.exon
    df.exon.all$hl = hl.matched[match(rownames(df.exon.all),names(hl.matched))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    df.exon.all$m = s.m[rownames(df.exon.all)]
    df.exon.all = na.omit(df.exon.all)
    df = df.exon.all
    
    HH = which(df$m>quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    HL = which(df$m>quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    LH = which(df$m<quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    LL = which(df$m<quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    Hx = c(HH,HL)
    Lx = c(LH,LL)
    xH = c(HH,LH)
    xL = c(LL,HL)
    boxplot(list(a = df$top_precision[Hx],b = df$top_precision[Lx]),
            names = c("High","Low"),ylab = "Top precision",main = paste(sample_name, "Exon: exp",sep = " "))
    boxplot(list(a = df$top_precision[xH],b = df$top_precision[xL]),
            names = c("High","Low"),ylab = "Top precision",main = paste(sample_name, "Exon: mhl",sep = " "))
    
    print(paste(sample_name,":",sep = ""))
    p1 = wilcox.test(df$top_precision[c(HH,HL)],df$top_precision[c(LH,LL)],alternative = "greater")$p.value
    print(paste("Exon   exp, pvalue:",format(p1,scientific = T,digits = 3),sep = ""))
    p2 = wilcox.test(df$top_precision[c(HH,LH)],df$top_precision[c(HL,LL)],alternative = "less")$p.value
    print(paste("Exon   mhl, pvalue:",format(p2,scientific = T,digits = 3),sep = ""))
    
    # mRNA half-life - better proportion regression
    df.inex.all = df.inex
    df.inex.all$hl = hl.matched[match(rownames(df.inex.all),names(hl.matched))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    df.inex.all$m = u.m[rownames(df.inex.all)]
    df.inex.all = na.omit(df.inex.all)
    df = df.inex.all
    
    HH = which(df$m>quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    HL = which(df$m>quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    LH = which(df$m<quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    LL = which(df$m<quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    Hx = c(HH,HL)
    Lx = c(LH,LL)
    xH = c(HH,LH)
    xL = c(LL,HL)
    boxplot(list(a = df$top_precision[Hx],b = df$top_precision[Lx]),
            names = c("High","Low"),ylab = "Top precision",main = paste(sample_name, "Intron: exp",sep = " "))
    boxplot(list(a = df$top_precision[xH],b = df$top_precision[xL]),
            names = c("High","Low"),ylab = "Top precision",main = paste(sample_name, "Intron: mhl",sep = " "))
    
    p3 = wilcox.test(df$top_precision[c(HH,HL)],df$top_precision[c(LH,LL)],alternative = "greater")$p.value
    print(paste("Intron exp, pvalue:",format(p3,scientific = T,digits = 3),sep = ""))
    p4 = wilcox.test(df$top_precision[c(HH,LH)],df$top_precision[c(HL,LL)],alternative = "less")$p.value
    print(paste("Intron mhl, pvalue:",format(p4,scientific = T,digits = 3),sep = ""))
    
    df.pval = rbind(df.pval,c(p1,p2,p3,p4))
    # summary for pvalue
    if(F){
      rownames(df.pval) = sample_names[sample_ids]
      # rename
      rownames(df.pval) = c("hForebrain","mForebrain","BoneMarrow","endocrinogenesis")
      colnames(df.pval) = c("ex_exp","ex_mhl","in_exp","in_mhl")
      
      # barplot -- expression level
      df.plot = -log10(df.pval)
      barplot(t(as.matrix(df.plot[,c(1,3)])),beside = T,col = c("blue","red"),ylab = "-log10(pvalue)",
              main = "effect of expression level")
      abline(h = -log10(0.05),lty = 2)
      barplot(t(as.matrix(df.plot[,c(2,4)])),beside = T,col = c("blue","red"),ylab = "-log10(pvalue)",
              main = "effect of mRNA half-life")
      abline(h = -log10(0.05),lty = 2)
      
      # barplot(t(as.matrix(df.plot[c(4,3,2,1),c(1,3)])),beside = T,col = c("blue","red"),xlab = "-log10(pvalue)",
      #         main = "effect of expression level",horiz = T)
      # barplot(t(as.matrix(df.plot[c(4,3,2,1),c(2,4)])),beside = T,col = c("blue","red"),xlab = "-log10(pvalue)",
      #         main = "effect of mRNA half-life",horiz = T)

    }
  }
  # investigate the effect of mean count and mRNA half-life: random precision
  if(T){
    par(mfrow = c(1,4))
    # match gene names
    rownames(df.inex) = toupper(rownames(df.inex))
    rownames(df.exon) = toupper(rownames(df.exon))
    rownames(gt) = toupper(rownames(gt))
    colnames(gt) = toupper(colnames(gt))
    names(s.m) = toupper(names(s.m))
    names(u.m) = toupper(names(u.m))
    
    # mRNA half-life - better proportion regression
    df.exon.all = df.exon
    df.exon.all$hl = hl.matched[match(rownames(df.exon.all),names(hl.matched))]
    df.exon.all$diff = df.exon.all$AUPR - df.exon.all$random_precision
    df.exon.all$m = s.m[rownames(df.exon.all)]
    df.exon.all = na.omit(df.exon.all)
    df = df.exon.all
    
    HH = which(df$m>quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    HL = which(df$m>quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    LH = which(df$m<quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    LL = which(df$m<quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    Hx = c(HH,HL)
    Lx = c(LH,LL)
    xH = c(HH,LH)
    xL = c(LL,HL)
    boxplot(list(a = df$random_precision[Hx],b = df$random_precision[Lx]),
            names = c("High","Low"),ylab = "random_precision",main = paste(sample_name, "Exon: exp",sep = " "))
    boxplot(list(a = df$random_precision[xH],b = df$random_precision[xL]),
            names = c("High","Low"),ylab = "random_precision",main = paste(sample_name, "Exon: mhl",sep = " "))
    
    print(paste(sample_name,":",sep = ""))
    p1 = wilcox.test(df$random_precision[c(HH,HL)],df$random_precision[c(LH,LL)],alternative = "greater")$p.value
    print(paste("Exon   exp, pvalue:",format(p1,scientific = T,digits = 3),sep = ""))
    p2 = wilcox.test(df$random_precision[c(HH,LH)],df$random_precision[c(HL,LL)],alternative = "less")$p.value
    print(paste("Exon   mhl, pvalue:",format(p2,scientific = T,digits = 3),sep = ""))
    
    # mRNA half-life - better proportion regression
    df.inex.all = df.inex
    df.inex.all$hl = hl.matched[match(rownames(df.inex.all),names(hl.matched))]
    df.inex.all$diff = df.inex.all$AUPR - df.inex.all$random_precision
    df.inex.all$m = u.m[rownames(df.inex.all)]
    df.inex.all = na.omit(df.inex.all)
    df = df.inex.all
    
    HH = which(df$m>quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    HL = which(df$m>quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    LH = which(df$m<quantile(df$m,0.5)&df$hl>quantile(df$hl,0.5))
    LL = which(df$m<quantile(df$m,0.5)&df$hl<quantile(df$hl,0.5))
    Hx = c(HH,HL)
    Lx = c(LH,LL)
    xH = c(HH,LH)
    xL = c(LL,HL)
    boxplot(list(a = df$random_precision[Hx],b = df$random_precision[Lx]),
            names = c("High","Low"),ylab = "random_precision",main = paste(sample_name, "Intron: exp",sep = " "))
    boxplot(list(a = df$random_precision[xH],b = df$random_precision[xL]),
            names = c("High","Low"),ylab = "random_precision",main = paste(sample_name, "Intron: mhl",sep = " "))
    
    p3 = wilcox.test(df$random_precision[c(HH,HL)],df$random_precision[c(LH,LL)],alternative = "greater")$p.value
    print(paste("Intron exp, pvalue:",format(p3,scientific = T,digits = 3),sep = ""))
    p4 = wilcox.test(df$random_precision[c(HH,LH)],df$random_precision[c(HL,LL)],alternative = "less")$p.value
    print(paste("Intron mhl, pvalue:",format(p4,scientific = T,digits = 3),sep = ""))
    
    df.pval = rbind(df.pval,c(p1,p2,p3,p4))
    # summary for pvalue
    if(F){
      rownames(df.pval) = sample_names[sample_ids]
      # rename
      rownames(df.pval) = c("hForebrain","mForebrain","BoneMarrow","endocrinogenesis")
      colnames(df.pval) = c("ex_exp","ex_mhl","in_exp","in_mhl")
      
      # barplot -- expression level
      df.plot = -log10(df.pval)
      barplot(t(as.matrix(df.plot[,c(1,3)])),beside = T,col = c("blue","red"),ylab = "-log10(pvalue)",
              main = "effect of expression level")
      abline(h = -log10(0.05),lty = 2)
      barplot(t(as.matrix(df.plot[,c(2,4)])),beside = T,col = c("blue","red"),ylab = "-log10(pvalue)",
              main = "effect of mRNA half-life")
      abline(h = -log10(0.05),lty = 2)
      
      # barplot(t(as.matrix(df.plot[c(4,3,2,1),c(1,3)])),beside = T,col = c("blue","red"),xlab = "-log10(pvalue)",
      #         main = "effect of expression level",horiz = T)
      # barplot(t(as.matrix(df.plot[c(4,3,2,1),c(2,4)])),beside = T,col = c("blue","red"),xlab = "-log10(pvalue)",
      #         main = "effect of mRNA half-life",horiz = T)
      
    }
  }
  # investigate the effect of mean count and mRNA half-life: compare intron with exon
  if(F){
    par(mfrow = c(1,1))
    # match gene names
    rownames(df.inex) = toupper(rownames(df.inex))
    rownames(df.exon) = toupper(rownames(df.exon))
    rownames(gt) = toupper(rownames(gt))
    colnames(gt) = toupper(colnames(gt))
    names(s.m) = toupper(names(s.m))
    names(u.m) = toupper(names(u.m))
    
    # mRNA half-life - better proportion regression
    df.diff = df.exon
    df.diff$diff = df.inex$AUPR - df.exon$AUPR
    df.diff$tpdiff = df.inex$top_precision - df.exon$top_precision
    df.diff$tp.in = df.inex$top_precision
    df.diff$tp.ex = df.exon$top_precision
    
    # df.diff$rocdiff = df.inex$AUROC - df.exon$AUROC
    # df.diff$roc.in = df.inex$AUROC
    # df.diff$roc.ex = df.exon$AUROC
    # id = which(df.diff$roc.ex<0.5&df.diff$roc.in<0.5)
    
    # remove useless target genes: both AUROC < 0.5
    # id = which(df.diff$tp.in<df.diff$random_precision&df.diff$tp.ex<df.diff$random_precision)
    # length(id)
    # df.diff = df.diff[-id,]
    
    df.diff$hl = hl.matched[match(rownames(df.diff),names(hl.matched))]
    df.diff$sm = s.m[rownames(df.diff)]
    df.diff$um = u.m[rownames(df.diff)]
    df.diff = na.omit(df.diff)
    df = df.diff
    
    # correlation between 3 factors
    cor(cbind(log2(df$hl),log2(df$sm),log2(df$um)))
    # plot(log2(df$sm),log2(df$um))
    
    
    # calculate proportion
    cal_prop = function(df,x){
      n_in = length(which(df$tpdiff[x]>0))
      n_ex = length(which(df$tpdiff[x]<0))
      prop = n_in/(n_in+n_ex)
      return(prop)
    }
    
    # effect of mRNA mean count
    H = which(df$sm>quantile(df$sm,0.5))
    L = which(df$sm<quantile(df$sm,0.5))
    prop.H = cal_prop(df,H)
    prop.L = cal_prop(df,L)
    barplot(c(prop.H,prop.L),names.arg = c("H","L"),main = "Effect of mRNA expression",ylim = c(0,1),
            ylab = "proportion: intron>exon")
    # boxplot(list(a = df$diff[H],b = df$diff[L]))
    p1 = wilcox.test(df$tpdiff[H],df$tpdiff[L],alternative = "greater")$p.value
    
    # effect of mRNA half-life
    H = which(df$hl>quantile(df$hl,0.5))
    L = which(df$hl<quantile(df$hl,0.5))
    prop.H = cal_prop(df,H)
    prop.L = cal_prop(df,L)
    barplot(c(prop.H,prop.L),names.arg = c("H","L"),main = "Effect of mRNA half-life",ylim = c(0,1),
            ylab = "proportion: intron>exon")
    # boxplot(list(a = df$diff[H],b = df$diff[L]))
    p2 = wilcox.test(df$tpdiff[H],df$tpdiff[L],alternative = "less")$p.value
    
    # effect of premRNA expression
    H = which(df$um>quantile(df$um,0.5))
    L = which(df$um<quantile(df$um,0.5))
    prop.H = cal_prop(df,H)
    prop.L = cal_prop(df,L)
    barplot(c(prop.H,prop.L),names.arg = c("H","L"),main = "Effect of pre-mRNA expression",ylim = c(0,1),
            ylab = "proportion: intron>exon")
    # boxplot(list(a = df$diff[H],b = df$diff[L]))
    p3 = wilcox.test(df$tpdiff[H],df$tpdiff[L],alternative = "greater")$p.value
    
    # # effect of premRNA expression
    # thrs = seq(0.1,0.9,0.1)
    # props.H = props.L = NULL
    # for(thr in thrs){
    #   H = which(df$um>quantile(df$um,thr))
    #   L = which(df$um<quantile(df$um,thr))
    #   prop.H = cal_prop(df,H)
    #   prop.L = cal_prop(df,L)
    # 
    #   props.H = c(props.H,prop.H)
    #   props.L = c(props.L,prop.L)
    # }
    # plot(thrs,props.H,type = "b",ylim = c(0,1),col = "green",ylab = "proportion: intron>exon")
    # points(thrs,props.L,type = "b",ylim = c(0,1),col = "grey")
    # abline(h = 0.5,col ="black",lty = 2)
    
    df.pval = rbind(df.pval,c(p1,p2,p3))
    # summary for pvalue
    if(F){
      rownames(df.pval) = sample_names[sample_ids]
      # rename
      rownames(df.pval) = c("hForebrain","mForebrain","BoneMarrow","endocrinogenesis")
      colnames(df.pval) = c("exp.s","mhl","exp.u")
      
      # barplot -- expression level
      df.plot = -log10(df.pval)
      barplot(t(as.matrix(df.plot)),beside = T,ylab = "-log10(pvalue)",
              main = "effect of 3 factors")
      abline(h = -log10(0.05),lty = 2)
      
    }
  }
}

}

dev.off()
