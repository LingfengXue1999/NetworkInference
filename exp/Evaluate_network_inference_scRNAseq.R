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
#Main function: calculate the results of A1 with threshold, and compare to Agr1 to calculate precision, recall and FPR
Evaluation <- function(A1,Agr1,thr){    #calculate the results of A1 with threshold, and compare to Agr1 to calculate precision, recall and FPR
  B = matrix(0,nrow=nrow(A1),ncol = ncol(A1)) 
  for(i in 1:nrow(A1)){
    for(j in 1:ncol(A1)){
      if(A1[i,j]>=thr) B[i,j]<- 1
    }
  } 
  a=b=c=d=0 
  for(i in 1:nrow(B)){
    for(j in 1:ncol(B)){
      if(i!=j){				#attention that here ignored the autoregulation
        if(Agr1[i,j]==1&&B[i,j]==1){a=a+1
        }else if(Agr1[i,j]==1&&B[i,j]==0) { b=b+1
        }else if(Agr1[i,j]==0&&B[i,j]==1) { c=c+1
        }else if(Agr1[i,j]==0&&B[i,j]==0) { d=d+1}
      }
    }
  }	
  if((a+c)==0){
    return(c(0,0,c/(c+d)))
  }else{
    return(c(a/(a+c),a/(a+b),c/(c+d))) ##c(precision,recall,FPR)
  }
}

# Evaluate with motif or dorothea
gt_type_names =c("dorothea","motif")
for(gt_type in 1:2){
  # dataframe contains: sample name, ground-truth name, input, EPR, AUROC, AUPR, random precision
  df = NULL
  data.pr = list()
  pr.id = 1
  for(sample_id in 1:10){
    # load NI result 
    # load network
    if(T){
      sample_names = c("hBCell","hCortOrga","hLiverPro","hLiverHomo","hNKCell","hCD4TCell","hMonoCell","hPlacenta","hSpleen","hTCell")
      sample_sources = rep("human",10)
      sample_name = sample_names[sample_id]
      sample_source = sample_sources[sample_id]
      
      # load exon result
      load(file = paste("./more_datasets_0301/results/genie3_exon_",sample_name,".RData",sep =""))
      exon.mat = weightMatrix
      # load inex result
      load(file = paste("./more_datasets_0301/results/genie3_inex_",sample_name,".RData",sep =""))
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
    
    # gt <--> sample
    if(gt_type ==1){
      if(sample_source=="human"){
        gt_id = 1
      }
      if(sample_source=="mouse"){
        gt_id = 4
      }
    }
    if(gt_type ==2){
      if(sample_source=="human"){
        gt_id = 6
      }
      if(sample_source=="mouse"){
        gt_id = 7
      }
    }
    
    # ground-truth network
    if(T){
      gt_names = c("dorothea","RegM","TRRUST_M","Doro_M","motif_sel","motif","motif_m")
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
    }
    
    # Evaluation
    print(paste(sample_name,gt_name,"Start!!!",sep = " "))
    mat = exon.mat[TFs,targets]
    if(T){
      A1 = mat
      #construct Agr matrix, regulation exists:1, non-exists:0
      Agr1 <- gt
      
      #calculate precisino, recall and FPR
      interval= 100
      precision = recall = FPR =  rep(0,times = interval+1)
      
      for(i in 1:(interval+1)){			#equal interval for threshold
        thr = (i-1 )/interval*max(A1)
        result = Evaluation(A1,Agr1,thr)
        precision[i] = result[1]
        recall[i] = result[2]
        FPR[i] = result[3]
        print(i/interval)
      }
      w <- cbind((0:interval)/interval*max(A1),precision,recall)		#compare precision and recall in different threshold
      recall.exon = recall
      precision.exon = precision
      FPR.exon = FPR
      
      # plot(recall,precision,type = 'b',main = paste(sample_name,gt_name,"exon",sep = " "))  #, ylim = c(0,0.65)
      # abline(h = precision[1], col = 'blue')
      #plot(FPR,recall)
      #calculate AUPR and AUROC with rectangle area
      
      # plot with number
      if(F){
        plot(recall,precision,type = 'b',col = "red",
             ylim = c(0,max(precision,precision)), xlim = c(0,0.1))
        title(main = paste(sample_name,gt_name,sep = " "),line = 2.5)
        abline(h = precision[1], col = 'blue')
        par(new = T)
        gt_T_num  =sum(as.vector(gt))
        plot(recall*gt_T_num,precision,type = 'b',
             col = "red", ylim = c(0,max(precision,precision)), xlim = c(0,0.1*gt_T_num),
             xaxt = "n", yaxt = "n", ylab = "", xlab = "")
        axis(side = 3)
      }
      AUPR = 0
      for(i in 1:(interval)){AUPR = AUPR+ (recall[i] - recall[i+1])*(precision[i]+precision[i+1])/2}
      AUROC = 0
      for(i in 1:interval){AUROC = AUROC+ (FPR[i] - FPR[i+1])*(recall[i]+recall[i+1])/2}
      random_precision = precision[1]
      
      # calculate EPR
      density_gt = sum(gt)/length(gt)
      thr = quantile(as.vector(A1),1-density_gt)
      result = Evaluation(A1,Agr1,thr)
      EP = result[1]   #precision
      EPR = EP / precision[1]
      # points(result[1],result[2],col = 'red',pch =16,cex = 2)
      x.exon = result[1]
      y.exon = result[2]
      
      # df:  sample name, ground-truth name, input, EPR, AUROC, AUPR, EP, random precision
      df.tmp = c(sample_name,gt_name,"exon", EPR, AUROC, AUPR, EP, random_precision)
      df = rbind(df,df.tmp)
      
      print(paste(sample_name,gt_name,"exon","Done!!!",sep = " "))
    }
    
    
    # Evaluation for inex.mat
    mat = inex.mat[TFs,targets]
    if(T){
      A1 = mat
      #construct Agr matrix, regulation exists:1, non-exists:0
      Agr1 <- gt
      
      #calculate precisino, recall and FPR
      interval= 100
      precision = recall = FPR =  rep(0,times = interval+1)
      # thr =  (0:interval )/interval*max(A1)
      # result = sapply(thr,Evaluation,A1 = A1, Agr1 = Agr1)
      
      for(i in 1:(interval+1)){			#equal interval for threshold
        thr = (i-1 )/interval*max(A1)
        result = Evaluation(A1,Agr1,thr)
        precision[i] = result[1]
        recall[i] = result[2]
        FPR[i] = result[3]
        print(i/interval)
      }
      w <- cbind((0:interval)/interval*max(A1),precision,recall)		#compare precision and recall in different threshold
      recall.inex = recall
      precision.inex = precision
      FPR.inex = FPR
      
      # plot(recall.inex,precision.inex,type = 'b',main = paste(sample_name,gt_name,sep = " "),col = "red",
      #      ylim = c(0,max(precision.inex,precision.exon)))
      # abline(h = precision[1], col = 'blue')
      # points(recall.exon,precision.exon,type = "b",col = "black")
      # # abline(h = precision.exon[1], col = 'blue')
      # 
      # 
      # plot(recall.inex,precision.inex,type = 'b',main = paste(sample_name,gt_name,sep = " "),col = "red",
      #      ylim = c(0,max(precision.inex,precision.exon)), xlim = c(0,0.1))
      # abline(h = precision[1], col = 'blue')
      # points(recall.exon,precision.exon,type = "b",col = "black")
      
      
      # plot with number
      if(F){
        plot(recall.inex,precision.inex,type = 'b',col = "red",log = "x",
             ylim = c(0,max(precision.inex,precision.exon)), xlim = c(1e-6,1),xlab = "recall",ylab = "precision")
        # title(main = paste(sample_name,gt_name,sep = " "),line = 2.5)
        title(paste(sample_name," PR curve", sep = ""))
        abline(h = precision[1], col = 'grey',lty = 2)
        points(recall.exon,precision.exon,type = "b",col = "black")
        # par(new = T)
        # gt_T_num  =sum(as.vector(gt))
        # plot(recall.inex*gt_T_num,precision.inex,type = 'b',
        #      col = "red", ylim = c(0,max(precision.inex,precision.exon)), xlim = c(0,0.1*gt_T_num),
        #      xaxt = "n", yaxt = "n", ylab = "", xlab = "")
        # axis(side = 3)
        cols = c("red","black","grey")
        ltys = c(1,1,2)
        legend("topright",c("inex","exon","random"),text.col = cols,col = cols,lty = ltys,inset = 0.05,cex = 1.3)
        
      }
      
      
      # plot(recall,precision,type = 'b',main = paste(sample_name,gt_name,"inex",sep = " "))  #, ylim = c(0,0.65)
      # abline(h = precision[1], col = 'blue')
      # plot(recall,precision,type = 'b',xlim = c(0,0.1),main = paste(sample_name,gt_name,"inex",sep = " "))  #, ylim = c(0,0.65)
      # abline(h = precision[1], col = 'blue')
      #plot(FPR,recall)
      #calculate AUPR and AUROC with rectangle area
      AUPR = 0
      for(i in 1:(interval)){AUPR = AUPR+ (recall[i] - recall[i+1])*(precision[i]+precision[i+1])/2}
      AUROC = 0
      for(i in 1:interval){AUROC = AUROC+ (FPR[i] - FPR[i+1])*(recall[i]+recall[i+1])/2}
      random_precision = precision[1]
      
      # calculate EPR
      density_gt = sum(gt)/length(gt)
      thr = quantile(as.vector(A1),1-density_gt)
      result = Evaluation(A1,Agr1,thr)
      EP = result[1]   #precision
      EPR = EP / precision[1]
      x.inex = result[1]
      y.inex = result[2]
      # points(x.inex,y.inex,col = 'red',pch =16,cex = 2)
      # points(x.exon,y.exon,col = 'black',pch =16,cex = 2)
      
      
      # df:  sample name, ground-truth name, input, EPR, AUROC, AUPR, EP, random precision
      df.tmp = c(sample_name,gt_name,"inex", EPR, AUROC, AUPR, EP, random_precision)
      df = rbind(df,df.tmp)
      
      print(paste(sample_name,gt_name,"inex","Done!!!",sep = " "))
      
    }
    
    data.pr[[pr.id]] = data.frame(recall_inex = recall.inex, recall_exon = recall.exon,
                                  precision_inex = precision.inex, precision_exon = precision.exon,
                                  FPR_inex = FPR.inex, FPR_exon = FPR.exon)
    names(data.pr)[pr.id] = sample_name
    pr.id = pr.id + 1
  }
  df = as.data.frame(df)
  colnames(df) = c("sample name", "ground-truth name", "input", "EPR", "AUROC", "AUPR", "EP", "random precision")
  rownames(df) = NULL
  save(df,data.pr,file = paste("./Evaluation_result_",gt_type_names[gt_type],".RData",sep = ""))
  
  
  
}


