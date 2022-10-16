# Evaluate network inference results for hESC data
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
# Evaluate for specific TF
Evaluation_TF = function(mat_targets,gt_targets){
  weights = sort(mat_targets)
  gt_targets = intersect(names(mat_targets),gt_targets)
  df = NULL
  for(j in 1:length(weights)){
    thr = weights[j]
    infer_targets = names(mat_targets)[which(mat_targets>=thr)]
    true_targets = intersect(gt_targets,infer_targets)
    TP = length(true_targets)
    FP = length(infer_targets) - length(true_targets)
    FN = length(gt_targets) - length(true_targets)
    TN = length(names(mat_targets)) - TP - FP - FN
    precision = TP/(TP+FP)
    recall = TP/(TP+FN)
    FPR = FP/(FP+TN)
    df = rbind(df,c(thr,precision,recall,FPR))
  }
  colnames(df) = c("thr","precision","recall","FPR")
  rownames(df) = NULL
  df = as.data.frame(df)
  
  return(df)
}
Evaluation_all_TF = function(mat,gt){
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
library(ggplot2)


# plot AUPR, AUROC, precision-recall curve, ROC for scRNA-seq datasets
gt_type_names =c("dorothea","motif")
for(gt_type in 1:2){
  load(file =  paste("./0114_Evaluation_result_",gt_type_names[gt_type],".RData",sep = ""))
  df.all = df
  data.pr.all = data.pr
  load(file =  paste("./more_datasets_0301/Evaluation_result_",gt_type_names[gt_type],".RData",sep = ""))
  df.all = rbind(df.all,df)
  data.pr.all = c(data.pr.all,data.pr)
  
  df = df.all
  data.pr = data.pr.all
  # save(df,data.pr,file = paste("./more_datasets_0301/Evaluation_results_30_samples/Evaluation_result_",gt_type_names[gt_type],".RData",sep = ""))
  
  colnames(df) = c("sample_name", "gt_ name", "input", "EPR", "AUROC", "AUPR", "EP", "random_precision")
  df = as.data.frame(as.matrix(df))
  df$EPR = as.numeric(as.character(df$EPR))
  df$AUROC = as.numeric(as.character(df$AUROC))
  df$AUPR = as.numeric(as.character(df$AUPR))
  df$EP = as.numeric(as.character(df$EP))
  df$random_precision = as.numeric(as.character(df$random_precision))
  
  # add species information
  species_info = c(rep("mouse",8),rep("human",6),rep("mouse",2),rep("human",14),rep("mouse",10),rep("human",20))
  df$species_info = species_info
  
  
  # change name of input
  dict = data.frame(inex = "intron",exon = "exon")
  df$input = factor(as.character(as.matrix(dict[as.character(df$input)])),levels = as.character(as.matrix(dict)) )
  # sort dataframe according to species
  df = df[order(df$species_info),]
  df$sample_name = factor(df$sample_name,levels = unique(df$sample_name))
  
  
  # plot AUROC
  p = ggplot() + geom_bar(data = df, aes(x = sample_name,y = AUROC,  fill = input),stat = "identity",
                          position = "dodge")+ scale_fill_manual(values=c("#F8766D","#00BFC4"))
  
  p = p + ylim(0,1) + labs(x = NULL,y = "AUROC",title = paste("Evaluate with ",gt_type_names[gt_type],sep = ""))
  p = p + theme(plot.title = element_text(size = 20)) + geom_hline(aes(yintercept=0.5),colour="blue", linetype="dashed")
  print(p)
  p1 = p
  
  # add random rows
  nrows = nrow(df)
  sub_rows = 2*(1:(nrows/2)) # choose even rows
  rows.add = df[sub_rows,]
  rows.add$AUPR = rows.add$random_precision
  rows.add$input = "random"
  df.add = rbind(df,rows.add)
  # plot AUPRC
  p = ggplot() + geom_bar(data = df.add, aes(x = sample_name,y = AUPR,  fill = input),stat = "identity",
                          position = "dodge")+ scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
  
  p = p + ylim(0,max(df.add$AUPR)) + labs(x = NULL,y = "AUPR",title = paste("Evaluate with ",gt_type_names[gt_type],sep = ""))
  p = p + theme(plot.title = element_text(size = 20))
  p2 = p
  
  #names(data.pr) = unique(df[,1])
  
  df.in = df[which(df$input=="intron"),]
  df.ex = df[which(df$input=="exon"),]
  
  
  
  
  boxplot(df.in$AUPR/df.ex$AUPR,log = "y")
  t.test(log2(df.in$AUPR/df.ex$AUPR))
  

  
  
  # plot detailed curves
  # plot precision-recall curve
  
  pdf(file = paste("./Evaluation_PRC_",gt_type_names[gt_type],".pdf",sep = ""), width = 10,height = 6)
  
  par(mfrow = c(2,4))
  plot_ids = 1:10
  for(i in plot_ids){
    sample_name = names(data.pr)[i]
    
    a = data.pr[[i]]
    recall.inex = a[,1]
    recall.exon = a[,2]
    precision.inex = a[,3]
    precision.exon = a[,4]
    if(gt_type ==1){xmin =1e-3}
    if(gt_type ==2){xmin =1e-6}
    plot(recall.inex,precision.inex,type = 'l',col = "red",log = "x",
         ylim = c(0,1), xlim = c(xmin,1),xlab = "recall",ylab = "precision")
    # title(main = paste(sample_name,gt_name,sep = " "),line = 2.5)
    title(paste(sample_name," PR curve",sep = ""))
    
    abline(h = precision.exon[1], col = 'grey',lty = 2)
    points(recall.inex,precision.inex,type = "l",col = "red",pch = 1)
    points(recall.exon,precision.exon,type = "l",col = "black",pch = 1)
    
  }
  dev.off()
  
  # plot ROC
  for(i in plot_ids){
    sample_name = names(data.pr)[i]
    
    a = data.pr[[i]]
    recall.inex = a[,1]
    recall.exon = a[,2]
    FPR.inex = a[,5]
    FPR.exon = a[,6]
    plot(FPR.inex,recall.inex,type = 'l',col = "red",
         ylim = c(0,1), xlim = c(0,1),xlab = "FPR",ylab = "TPR")
    
    # title(main = paste(sample_name,gt_name,sep = " "),line = 2.5)
    title(paste(sample_name," ROC curve",sep = ""))
    points(FPR.exon,recall.exon,type = "l",col = "black",pch = 1)
    
    abline(a = 0, b = 1, col = 'grey',lty = 2)
  }
  
  # plot top precision for motif evaluation
  if(gt_type==1){
    # calculate top precision: average precision at the limit when recall<0.01
    df.tp = NULL
    for(i in 1:length(data.pr)){
      a = data.pr[[i]]
      recall.inex = a[,1]
      recall.exon = a[,2]
      precision.inex = a[,3]
      precision.exon = a[,4]
      cal_tp = function(precision,recall,thr){
        id = which(recall<thr&recall>0)
        AUPR = 0
        for(i in id[-length(id)]){AUPR = AUPR+ (recall[i] - recall[i+1])*(precision[i]+precision[i+1])/2}
        AUPR = AUPR/(recall[id[1]]-recall[id[length(id)]])
        return(AUPR)
      }
      thr = 0.1
      tp.inex = cal_tp(precision.inex,recall.inex,thr)
      tp.exon = cal_tp(precision.exon,recall.exon,thr)
      
      df.tp = rbind(df.tp,c(tp.inex,tp.exon,precision.inex[1]))
    }
    df.tp = as.data.frame(df.tp)
    colnames(df.tp) = c("tp.inex","tp.exon","rand")
    rownames(df.tp) = names(data.pr)
    # sort df.tp
    df.tp = df.tp[as.vector(unique(df$sample_name)),]
    
    # plot results
    df.plot = data.frame(tp = c(df.tp[,1],df.tp[,2],df.tp[,3]),input = c(rep("intron",nrow(df.tp)),rep("exon",nrow(df.tp)),rep("random",nrow(df.tp))),
                         sample_name = rep(rownames(df.tp),3))
    df.plot$input = factor(df.plot$input,levels= c("intron","exon","random"))
    # sample_names = rownames(df.tp)[c(1,8,2,5,3,4,6,7)]
    # df.plot$sample_name = factor(df.plot$sample_name,levels = sample_names)
    df.plot$sample_name = factor(df.plot$sample_name,levels = unique(df.plot$sample_name))
    p = ggplot() + geom_bar(data = df.plot, aes(x = sample_name,y = tp,  fill = input),stat = "identity",
                            position = "dodge")+ scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
    
    p = p + ylim(0,0.5) + labs(x = NULL,y = "Average Early Precision",title = paste("Evaluate with ",gt_type_names[gt_type],sep = ""))
    p = p + theme(plot.title = element_text(size = 20))
    
  }
  if(gt_type==2){
    # calculate top precision: average precision at the limit when recall<0.01
    df.tp = NULL
    for(i in 1:length(data.pr)){
      a = data.pr[[i]]
      recall.inex = a[,1]
      recall.exon = a[,2]
      precision.inex = a[,3]
      precision.exon = a[,4]
      cal_tp = function(precision,recall,thr){
        id = which(recall<thr&recall>0)
        AUPR = 0
        for(i in id[-length(id)]){AUPR = AUPR+ (recall[i] - recall[i+1])*(precision[i]+precision[i+1])/2}
        AUPR = AUPR/(recall[id[1]]-recall[id[length(id)]])
        return(AUPR)
      }
      thr = 0.01
      tp.inex = cal_tp(precision.inex,recall.inex,thr)
      tp.exon = cal_tp(precision.exon,recall.exon,thr)
      
      df.tp = rbind(df.tp,c(tp.inex,tp.exon,precision.inex[1]))
    }
    df.tp = as.data.frame(df.tp)
    colnames(df.tp) = c("tp.inex","tp.exon","rand")
    rownames(df.tp) = names(data.pr)
    # sort df.tp
    df.tp = df.tp[as.vector(unique(df$sample_name)),]
    
    # plot results
    df.plot = data.frame(tp = c(df.tp[,1],df.tp[,2],df.tp[,3]),input = c(rep("intron",nrow(df.tp)),rep("exon",nrow(df.tp)),rep("random",nrow(df.tp))),
                         sample_name = rep(rownames(df.tp),3))
    df.plot$input = factor(df.plot$input,levels= c("intron","exon","random"))
    df.plot$sample_name = factor(df.plot$sample_name,levels = unique(df.plot$sample_name))
    p = ggplot() + geom_bar(data = df.plot, aes(x = sample_name,y = tp,  fill = input),stat = "identity",
                            position = "dodge")+ scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
    
    p = p + ylim(0,0.4) + labs(x = NULL,y = "Average Early Precision",title = paste("Evaluate with ",gt_type_names[gt_type],sep = ""))
    p = p + theme(plot.title = element_text(size = 20))
    print(p)
  }
  p3 = p
  dev.off()
  
  # plot 3 results
  pdf(file = paste("./Evaluation_summary_",gt_type_names[gt_type],".pdf",sep = ""), width = 20,height = 4)
  
  print(p2)
  print(p1)
  print(p3)
  dev.off()
  
  
  pdf(file = paste("./Evaluation_summary2_",gt_type_names[gt_type],".pdf",sep = ""), width = 10,height = 5)
  par(mfrow = c(1,2))
  plot(df.ex$AUPR,df.in$AUPR,xlab = "AUPR for exon method",ylab = "AUPR for intron method", main = "Compare intron with exon",
       xlim = c(0,0.13),ylim = c(0,0.13))
  abline(a = 0, b = 1,col = "green")
  
  plot(df.ex$AUROC,df.in$AUROC,xlab = "AUROC for exon method",ylab = "AUROC for intron method", main = "Compare intron with exon",
       xlim = c(0.4,0.8),ylim = c(0.4,0.8))
  abline(a = 0, b = 1,col = "green")
  
  plot(df.tp$tp.exon,df.tp$tp.inex,xlab = "Average early precision for exon method",ylab = "Average early precision for intron method", main = "Compare intron with exon",
       xlim = c(0,0.5),ylim = c(0,0.5))
  abline(a = 0, b = 1,col = "green")
  dev.off()
  
  
  
  
}


# Evaluate each TF
pdf(file = "./Evaluate_each_TF_hFB.pdf",width = 8, height = 4.5)
par(mfrow = c(1,2))
TFs = c("NEUROD2","SOX2","SOX2")
for(TF_id in 1:3){
  TF = TFs[TF_id]
  gt_id = TF_id
  
  # load inferred network
  sample_id = 1
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
  # compare intron and exon -> change size
  if(T){
    dim(exon.mat)
    dim(inex.mat)
    overlap_targets = intersect(colnames(exon.mat),colnames(inex.mat))
    exon.mat = exon.mat[,overlap_targets]
    inex.mat = inex.mat[,overlap_targets]
  }
  # load ground-truth
  if(T){
    # NEUROD2 ChIP targets
    if(gt_id ==1){
      gt = read.csv(file = "./Extract_network_links/ChIP_data/NEUROD2_targets.csv",header = T)
      top_num = 1045
      gt = gt[1:top_num,]
      gt_targets = toupper(gt[,1])
      gt_name = "ChIP"
    }
    # Dorothea data
    if(gt_id ==2){
      load(file ="./dorothea/dorothea_database.RData")
      tmp = unique(dorothea_hs_high_confidence)
      tmp = as.matrix(tmp)
      reference = cbind(as.character(tmp[,1]),as.character(tmp[,3]))
      gt = list_to_mat(reference)
      
      gt_targets = toupper(colnames(gt)[which(gt[TF,]==1)])
      gt_name = "Dorothea"
    }
    # motif network
    if(gt_id==3){
      load(file = "./Rcistarget/motif_link_and_mat.RData")
      gt = motif_mat
      
      gt_targets = toupper(colnames(gt)[which(gt[TF,]==1)])
      gt_name = "Motif"
    }
  }
  # plot precision-recall curve for the TF
  if(T){
    df = cbind(Evaluation_TF(inex.mat[TF,],gt_targets),Evaluation_TF(exon.mat[TF,],gt_targets))
    colnames(df) = c("thr.in","pre.in","rec.in","FPR.in","thr.ex","pre.ex","rec.ex","FPR.ex")
    plot(df$rec.in,df$pre.in,ylim = c(0,1),type = "b",col = "red",main = paste(TF,gt_name,sep = " "),
         xlab = "recall",ylab = "precision")
    points(df$rec.ex,df$pre.ex,type = "b",col = "black")
    
    # plot ROC for TF
    if(T){
      plot(df$FPR.in,df$rec.in,ylim = c(0,1),type = "b",col = "red",main = paste(TF,gt_name,sep = " "),
           xlab = "FPR",ylab = "TPR")
      points(df$FPR.ex,df$rec.ex,type = "b",col = "black")
      abline(a = 0,b = 1,col = "grey")
      
    }
  }
}
dev.off()

# plot results for all TFs
pdf(file = "./Evaluate_all_TFs_hFB_doro.pdf",width = 16, height = 9)
if(T){
  # load inferred network
  sample_id = 1
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
  # compare intron and exon -> change size
  if(T){
    dim(exon.mat)
    dim(inex.mat)
    overlap_targets = intersect(colnames(exon.mat),colnames(inex.mat))
    exon.mat = exon.mat[,overlap_targets]
    inex.mat = inex.mat[,overlap_targets]
  }
  # load dorothea
  if(T){
    load(file ="./dorothea/dorothea_database.RData")
    tmp = unique(dorothea_hs_high_confidence)
    tmp = as.matrix(tmp)
    reference = cbind(as.character(tmp[,1]),as.character(tmp[,3]))
    gt = list_to_mat(reference)
    
    gt_targets = toupper(colnames(gt)[which(gt[TF,]==1)])
    gt_name = "Dorothea"
  }
  # change size
  if(T){
    TFs = intersect(rownames(gt),rownames(inex.mat))
    targets = intersect(colnames(gt),colnames(inex.mat))
    gt = gt[TFs,targets]
    exon.mat = exon.mat[TFs,targets]
    inex.mat = inex.mat[TFs,targets]
  }
  
  df.exon = Evaluation_all_TF(exon.mat,gt)
  df.inex = Evaluation_all_TF(inex.mat,gt)
  
  # AUPR dataframe
  df.comp = data.frame(TF = rep(rownames(df.exon),3), AUPR = c(df.inex$AUPR,df.exon$AUPR,df.inex$random_precision),
                       data.type = c(rep("intron",nrow(df.exon)),rep("exon",nrow(df.exon)),rep("random",nrow(df.exon))),
                       n = rep(df.exon$n_true_targets,3))
  df.comp$data.type = factor(df.comp$data.type,levels = c("intron","exon","random"))
  df.comp$TF = factor(df.comp$TF, levels = unique(df.comp$TF[order(df.comp$n,decreasing = T)]))
  df.plot = df.comp[which(df.comp$n>=5),]
  p = ggplot(df.plot, aes(x = TF, y = AUPR, fill = data.type)) + labs(title = "Evaluate each TF",x = NULL)
  p = p +  geom_bar(stat = "identity", position =position_dodge())
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
  p = p + geom_text(mapping=aes(x=TF,y=0,label=paste(n,"\n","targets",sep="")),colour="black",position = position_dodge(0),size=5)
  print(p)
}
dev.off()
pdf(file = "./Evaluate_all_TFs_hFB_motif.pdf",width = 16, height = 9)
if(T){
  # load inferred network
  sample_id = 1
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
  # compare intron and exon -> change size
  if(T){
    dim(exon.mat)
    dim(inex.mat)
    overlap_targets = intersect(colnames(exon.mat),colnames(inex.mat))
    exon.mat = exon.mat[,overlap_targets]
    inex.mat = inex.mat[,overlap_targets]
  }
  # load Motif
  if(T){
    load(file = "./Rcistarget/motif_link_and_mat.RData")
    gt = motif_mat
  }
  # change size
  if(T){
    TFs = intersect(rownames(gt),rownames(inex.mat))
    targets = intersect(colnames(gt),colnames(inex.mat))
    gt = gt[TFs,targets]
    exon.mat = exon.mat[TFs,targets]
    inex.mat = inex.mat[TFs,targets]
  }
  
  # df.exon = Evaluation_all_TF(exon.mat,gt)
  # df.inex = Evaluation_all_TF(inex.mat,gt)
  # save(df.exon,df.inex,file = "./eval_all_TFs_hFB_motif.RData")
  load(file = "./eval_all_TFs_hFB_motif.RData")
  
  diff = df.inex$AUPR - df.exon$AUPR
  names(diff) = rownames(df.inex)
  diff = sort(diff,decreasing = T)
  n =10
  id.in = match(names(head(diff,n)),rownames(df.inex))
  id.ex = match(names(tail(diff,n)),rownames(df.inex))
  
  # AUPR: intron better
  id = id.in
  df.comp = data.frame(TF = rep(rownames(df.exon[id,]),3), AUPR = c(df.inex$AUPR[id],df.exon$AUPR[id],df.inex$random_precision[id]),
                       data.type = c(rep("intron",nrow(df.exon[id,])),rep("exon",nrow(df.exon[id,])),rep("random",nrow(df.exon[id,]))),
                       n = rep(df.exon$n_true_targets[id],3))
  df.comp$data.type = factor(df.comp$data.type,levels = c("intron","exon","random"))
  df.comp$TF = factor(df.comp$TF, levels = unique(df.comp$TF[order(df.comp$n,decreasing = T)]))
  df.plot = df.comp
  p = ggplot(df.plot, aes(x = TF, y = AUPR, fill = data.type)) + labs(title = "Evaluate each TF: intron better TFs",x = NULL)
  p = p +  geom_bar(stat = "identity", position =position_dodge())
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
  p = p + geom_text(mapping=aes(x=TF,y=0,label=paste(n,"\n","targets",sep="")),colour="black",position = position_dodge(0),size=5)
  print(p)
  # AUPR: exon better
  id = id.ex
  df.comp = data.frame(TF = rep(rownames(df.exon[id,]),3), AUPR = c(df.inex$AUPR[id],df.exon$AUPR[id],df.inex$random_precision[id]),
                       data.type = c(rep("intron",nrow(df.exon[id,])),rep("exon",nrow(df.exon[id,])),rep("random",nrow(df.exon[id,]))),
                       n = rep(df.exon$n_true_targets[id],3))
  df.comp$data.type = factor(df.comp$data.type,levels = c("intron","exon","random"))
  df.comp$TF = factor(df.comp$TF, levels = unique(df.comp$TF[order(df.comp$n,decreasing = T)]))
  df.plot = df.comp
  p = ggplot(df.plot, aes(x = TF, y = AUPR, fill = data.type)) + labs(title = "Evaluate each TF: exon better TFs",x = NULL)
  p = p +  geom_bar(stat = "identity", position =position_dodge())
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
  p = p + geom_text(mapping=aes(x=TF,y=0,label=paste(n,"\n","targets",sep="")),colour="black",position = position_dodge(0),size=5)
  print(p)
}
dev.off()
