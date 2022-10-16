# Network motif analysis for network inference # 2021/7/29
mat_to_link = function(a){
  TFs = rownames(a)
  targets = colnames(a)
  linklist = NULL
  for(i in 1:nrow(a)){
    for(j in 1:ncol(a)){
      linklist = rbind(linklist,c(TFs[i],targets[j],as.numeric(a[i,j])))
    }
  }
  linklist = data.frame(TFs = linklist[,1],targets = linklist[,2],weight = as.numeric(linklist[,3]))
  return(linklist)
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

Find_motifs = function(x){
  mat = as_adjacency_matrix(x)
  Find_fanin = function(x){
    # 2->1, 3->1
    df = NULL
    indegree = degree(x,mode = "in")
    id = which(indegree>=2)
    if(length(id)==0){
      return(NULL)
    }
    for(i in id){
      nodes = unlist(neighborhood(x,order = 1,nodes = V(x)[i],mode = "in",mindist = 1))
      # select results
      a = cbind(i,t(combn(nodes,2)))
      df = rbind(df,a)
    }
    if(is.null(df)){
      return(NULL)
    }else{
      link_num = apply(df,1,function(x){return(sum(sum(mat[x,x])))})
      df = df[which(link_num==2),]
    }
    colnames(df) = c("1","2","3")
    rownames(df) = NULL
    
    
    return(df)
  }
  Find_fanout = function(x){
    # 1->2, 1->3
    df = NULL
    outdegree = degree(x,mode = "out")
    id = which(outdegree>=2)
    if(length(id)==0){
      return(NULL)
    }
    for(i in id){
      nodes = unlist(neighborhood(x,order = 1,nodes = V(x)[i],mode = "out",mindist = 1))
      # select results
      a = cbind(i,t(combn(nodes,2)))
      df = rbind(df,a)
    }
    if(is.null(df)){
      return(NULL)
    }else{
      colnames(df) = c("1","2","3")
      link_num = apply(df,1,function(x){return(sum(sum(mat[x,x])))})
      df = df[which(link_num==2),]
    }
    rownames(df) = NULL
    
    
    return(df)
  }
  Find_cascade = function(x){
    # 1->2->3
    df = NULL
    outdegree = degree(x,mode = "out")
    id = which(outdegree>=1)
    if(length(id)==0){
      return(NULL)
    }
    for(i in id){
      nodes = unlist(neighborhood(x,order = 1,nodes = V(x)[i],mode = "out",mindist = 1))
      for(j in nodes){
        if(outdegree[j]<1){next}
        nodes3 = unlist(neighborhood(x,order = 1,nodes = V(x)[j],mode = "out",mindist = 1))
        a = cbind(i,j,nodes3)
        df = rbind(df,a)
      }
    }
    link_num = apply(df,1,function(x){return(sum(sum(mat[x,x])))})
    df = df[which(link_num==2),]
    if(is.null(df)){
      return(NULL)
    }
    colnames(df) = c("1","2","3")
    rownames(df) = NULL
    
    return(df)
  }
  Find_FFL = function(x){
    # 1->2, 1->3, 2->3
    df = NULL
    outdegree = degree(x,mode = "out")
    id = which(outdegree>=2)
    if(length(id)==0){
      return(NULL)
    }
    for(i in id){
      nodes = unlist(neighborhood(x,order = 1,nodes = V(x)[i],mode = "out",mindist = 1))
      for(j in nodes){
        nodes3 = unlist(neighborhood(x,order = 1,nodes = V(x)[j],mode = "out",mindist = 1))
        id3 = intersect(nodes3,nodes)
        if(length(id3)==0){next}
        a = cbind(i,j,id3)
        df = rbind(df,a)
      }
    }
    if(is.null(df)){
      return(NULL)
    }else{
      link_num = apply(df,1,function(x){return(sum(sum(mat[x,x])))})
      df = df[which(link_num==3),]
    }
    colnames(df) = c("1","2","3")
    rownames(df) = NULL
    return(df)
  }
  a = Find_fanin(x)
  b = Find_fanout(x)
  c = Find_cascade(x)
  d = Find_FFL(x)
  return(list(fanin=a,fanout=b,cascade=c,FFL=d))
}
library(igraph)
library(Matrix)
# load ground-truth network and inferred network
dir = 'dyngen/check_controllability/'

# plot Fan-out error
pdf("./motif_error_4backbones.pdf",width = 7, height = 7)
par(mfrow = c(2,2))
for(data_num in 1:4){
  net_id = 1
  rep_id = 0
  seed_id = 1
  if(T){
    load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
    
    net = model[["feature_network"]][,c(1,2,5)]
    # modify links: remain only TFs as regulators
    feature_info = model[["feature_info"]]
    TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
    net = net[which(is.element(net$from,TFs)),]
    # remove self regulation
    net = net[which(net$from!=net$to),]
    
    # load NI result 
    if(T){
      # load exon result
      load(file = paste(dir,'Network_inference/genie3_output_exon_',data_num,'_net',net_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
      exon.mat = weightMatrix
      
      # load inex result
      load(file = paste(dir,'Network_inference/genie3_output_inex_',data_num,'_net',net_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
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
    # mat to linklist
    inex.link = mat_to_link(inex.mat)
    inex.link = inex.link[order(inex.link$weight,decreasing = T),]
    inex.link$confidence = seq(1,0,length.out = nrow(inex.link))
    mat.conf.inex = list_to_mat(inex.link[,c(1,2,4)])
    
    exon.link = mat_to_link(exon.mat)
    exon.link = exon.link[order(exon.link$weight,decreasing = T),]
    exon.link$confidence = seq(1,0,length.out = nrow(exon.link))
    mat.conf.exon = list_to_mat(exon.link[,c(1,2,4)])
    
    # overlap genes
    TFs = rownames(inex.mat)
    targets = colnames(inex.mat)
    net = net[which(is.element(net$from,TFs)),]
    net = net[which(is.element(net$to,targets)),]
    
    
  }
  
  # identify all motifs
  x = graph_from_data_frame(d = net,directed = T)
  z = mat.conf.exon
  motifs = Find_motifs(x)
  print(data_num)
  print(unlist(lapply(motifs,nrow)))
  # print(motifs(x, 3)[c(3,7,5,8)]) 
  
  
  # determine prediction confidence of all instances for motifs
  # par(mfrow= c(2,2))
  motif_names = c("Fan_in","Fan_out","Cascade","FFL")
  input_names = c("inex","exon")
  
  for(motif_id in c(2)){
    motif_name = motif_names[motif_id]
    y = motifs[[motif_id]]
    for(input_id in c(1)){
      if(input_id == 1){
        z = mat.conf.inex
      }else{
        z = mat.conf.exon
      }
      input_name = input_names[input_id]
      
      # number of TFs in motif
      num_TFs = apply(y,1,function(y){return(length(intersect(names(V(x)[y]),TFs)))})
      table(num_TFs)
      TF_ids = match(TFs,names(V(x)))
      
      meta_list = NULL
      for(i in 1:3){
        for(j in 1:3){
          if(i!=j){
            meta_list = rbind(meta_list,c(i,j))
          }
        }
      }
      df.conf = NULL
      for(i in 1:nrow(y)){
        tmp = rep(NA,6)
        meta_sub = which(is.element(meta_list[,1],which(is.element(y[i,],TF_ids))))
        for(j in meta_sub){
          tmp[j] = z[names(V(x)[y[i,meta_list[j,1]]]), names(V(x)[y[i,meta_list[j,2]]])]
        }
        df.conf = rbind(df.conf,tmp)
      }
      df.conf = as.data.frame(df.conf)
      colnames(df.conf) = paste(meta_list[,1],meta_list[,2],sep = "")
      rownames(df.conf) = rownames(y)
      
      # determine the background prediction confidence
      true_egdes = net[,1:2]
      true_edges.conf = apply(net,1,function(x){return(z[x[1],x[2]])})
      absent_edges.conf = setdiff(as.numeric(z),true_edges.conf)
      
      # evaluating prediction confidence
      # C12, C13, C23
      if(motif_id==2){
        pvalue = wilcox.test(df.conf$`23`,absent_edges.conf)[["p.value"]]
        title = paste(data_num," ",input_name," ",motif_name," C23 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`23`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      if(motif_id==3){
        pvalue = wilcox.test(df.conf$`13`,absent_edges.conf)[["p.value"]]
        title = paste(data_num," ",input_name," ",motif_name," C13 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`13`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      boxplot(df.conf$`12`,df.conf$`13`,df.conf$`23`,true_edges.conf,absent_edges.conf,names = c("C12","C13","C23","True","False"),
              main = title)
    }
    tmp = c(data_num,motif_id,input_id,values)
    df.backbone = rbind(df.backbone,tmp)
  }
}
# plot Cascade error
for(data_num in 1:4){
  net_id = 1
  rep_id = 0
  seed_id = 1
  if(T){
    load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
    
    net = model[["feature_network"]][,c(1,2,5)]
    # modify links: remain only TFs as regulators
    feature_info = model[["feature_info"]]
    TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
    net = net[which(is.element(net$from,TFs)),]
    # remove self regulation
    net = net[which(net$from!=net$to),]
    
    # load NI result 
    if(T){
      # load exon result
      load(file = paste(dir,'Network_inference/genie3_output_exon_',data_num,'_net',net_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
      exon.mat = weightMatrix
      
      # load inex result
      load(file = paste(dir,'Network_inference/genie3_output_inex_',data_num,'_net',net_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
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
    # mat to linklist
    inex.link = mat_to_link(inex.mat)
    inex.link = inex.link[order(inex.link$weight,decreasing = T),]
    inex.link$confidence = seq(1,0,length.out = nrow(inex.link))
    mat.conf.inex = list_to_mat(inex.link[,c(1,2,4)])
    
    exon.link = mat_to_link(exon.mat)
    exon.link = exon.link[order(exon.link$weight,decreasing = T),]
    exon.link$confidence = seq(1,0,length.out = nrow(exon.link))
    mat.conf.exon = list_to_mat(exon.link[,c(1,2,4)])
    
    # overlap genes
    TFs = rownames(inex.mat)
    targets = colnames(inex.mat)
    net = net[which(is.element(net$from,TFs)),]
    net = net[which(is.element(net$to,targets)),]
    
    
  }
  
  # identify all motifs
  x = graph_from_data_frame(d = net,directed = T)
  z = mat.conf.exon
  motifs = Find_motifs(x)
  print(data_num)
  print(unlist(lapply(motifs,nrow)))
  # print(motifs(x, 3)[c(3,7,5,8)]) 
  
  
  # determine prediction confidence of all instances for motifs
  # par(mfrow= c(2,2))
  motif_names = c("Fan_in","Fan_out","Cascade","FFL")
  input_names = c("inex","exon")
  
  for(motif_id in c(3)){
    motif_name = motif_names[motif_id]
    y = motifs[[motif_id]]
    for(input_id in c(1)){
      if(input_id == 1){
        z = mat.conf.inex
      }else{
        z = mat.conf.exon
      }
      input_name = input_names[input_id]
      
      # number of TFs in motif
      num_TFs = apply(y,1,function(y){return(length(intersect(names(V(x)[y]),TFs)))})
      table(num_TFs)
      TF_ids = match(TFs,names(V(x)))
      
      meta_list = NULL
      for(i in 1:3){
        for(j in 1:3){
          if(i!=j){
            meta_list = rbind(meta_list,c(i,j))
          }
        }
      }
      df.conf = NULL
      for(i in 1:nrow(y)){
        tmp = rep(NA,6)
        meta_sub = which(is.element(meta_list[,1],which(is.element(y[i,],TF_ids))))
        for(j in meta_sub){
          tmp[j] = z[names(V(x)[y[i,meta_list[j,1]]]), names(V(x)[y[i,meta_list[j,2]]])]
        }
        df.conf = rbind(df.conf,tmp)
      }
      df.conf = as.data.frame(df.conf)
      colnames(df.conf) = paste(meta_list[,1],meta_list[,2],sep = "")
      rownames(df.conf) = rownames(y)
      
      # determine the background prediction confidence
      true_egdes = net[,1:2]
      true_edges.conf = apply(net,1,function(x){return(z[x[1],x[2]])})
      absent_edges.conf = setdiff(as.numeric(z),true_edges.conf)
      
      # evaluating prediction confidence
      # C12, C13, C23
      if(motif_id==2){
        pvalue = wilcox.test(df.conf$`23`,absent_edges.conf)[["p.value"]]
        title = paste(data_num," ",input_name," ",motif_name," C23 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`23`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      if(motif_id==3){
        pvalue = wilcox.test(df.conf$`13`,absent_edges.conf)[["p.value"]]
        title = paste(data_num," ",input_name," ",motif_name," C13 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`13`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      boxplot(df.conf$`12`,df.conf$`13`,df.conf$`23`,true_edges.conf,absent_edges.conf,names = c("C12","C13","C23","True","False"),
              main = title)
    }
    tmp = c(data_num,motif_id,input_id,values)
    df.backbone = rbind(df.backbone,tmp)
  }
}
df.backbone = NULL
for(data_num in 1:4){
  net_id = 1
  rep_id = 0
  seed_id = 1
  if(T){
    load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
    
    net = model[["feature_network"]][,c(1,2,5)]
    # modify links: remain only TFs as regulators
    feature_info = model[["feature_info"]]
    TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
    net = net[which(is.element(net$from,TFs)),]
    # remove self regulation
    net = net[which(net$from!=net$to),]
    
    # load NI result 
    if(T){
      # load exon result
      load(file = paste(dir,'Network_inference/genie3_output_exon_',data_num,'_net',net_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
      exon.mat = weightMatrix

      # load inex result
      load(file = paste(dir,'Network_inference/genie3_output_inex_',data_num,'_net',net_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
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
    # mat to linklist
    inex.link = mat_to_link(inex.mat)
    inex.link = inex.link[order(inex.link$weight,decreasing = T),]
    inex.link$confidence = seq(1,0,length.out = nrow(inex.link))
    mat.conf.inex = list_to_mat(inex.link[,c(1,2,4)])
    
    exon.link = mat_to_link(exon.mat)
    exon.link = exon.link[order(exon.link$weight,decreasing = T),]
    exon.link$confidence = seq(1,0,length.out = nrow(exon.link))
    mat.conf.exon = list_to_mat(exon.link[,c(1,2,4)])
    
    # overlap genes
    TFs = rownames(inex.mat)
    targets = colnames(inex.mat)
    net = net[which(is.element(net$from,TFs)),]
    net = net[which(is.element(net$to,targets)),]
    
    
  }
  
  # identify all motifs
  x = graph_from_data_frame(d = net,directed = T)
  z = mat.conf.exon
  motifs = Find_motifs(x)
  print(data_num)
  print(unlist(lapply(motifs,nrow)))
  # print(motifs(x, 3)[c(3,7,5,8)]) 
  
  
  # determine prediction confidence of all instances for motifs
  # par(mfrow= c(2,2))
  motif_names = c("Fan_in","Fan_out","Cascade","FFL")
  input_names = c("inex","exon")
  
  for(motif_id in c(2,3)){
    motif_name = motif_names[motif_id]
    y = motifs[[motif_id]]
    for(input_id in c(1)){
      if(input_id == 1){
        z = mat.conf.inex
      }else{
        z = mat.conf.exon
      }
      input_name = input_names[input_id]
      
      # number of TFs in motif
      num_TFs = apply(y,1,function(y){return(length(intersect(names(V(x)[y]),TFs)))})
      table(num_TFs)
      TF_ids = match(TFs,names(V(x)))
      
      meta_list = NULL
      for(i in 1:3){
        for(j in 1:3){
          if(i!=j){
            meta_list = rbind(meta_list,c(i,j))
          }
        }
      }
      df.conf = NULL
      for(i in 1:nrow(y)){
        tmp = rep(NA,6)
        meta_sub = which(is.element(meta_list[,1],which(is.element(y[i,],TF_ids))))
        for(j in meta_sub){
          tmp[j] = z[names(V(x)[y[i,meta_list[j,1]]]), names(V(x)[y[i,meta_list[j,2]]])]
        }
        df.conf = rbind(df.conf,tmp)
      }
      df.conf = as.data.frame(df.conf)
      colnames(df.conf) = paste(meta_list[,1],meta_list[,2],sep = "")
      rownames(df.conf) = rownames(y)
      
      # determine the background prediction confidence
      true_egdes = net[,1:2]
      true_edges.conf = apply(net,1,function(x){return(z[x[1],x[2]])})
      absent_edges.conf = setdiff(as.numeric(z),true_edges.conf)
      
      # evaluating prediction confidence
      # C12, C13, C23
      if(motif_id==2){
        pvalue = wilcox.test(df.conf$`23`,absent_edges.conf)[["p.value"]]
        title = paste(data_num," ",input_name," ",motif_name," C23 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`23`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      if(motif_id==3){
        pvalue = wilcox.test(df.conf$`13`,absent_edges.conf)[["p.value"]]
        title = paste(data_num," ",input_name," ",motif_name," C13 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`13`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      # boxplot(df.conf$`12`,df.conf$`13`,df.conf$`23`,true_edges.conf,absent_edges.conf,names = c("C12","C13","C23","True","False"),
      #         main = title)
    }
    tmp = c(data_num,motif_id,input_id,values)
    df.backbone = rbind(df.backbone,tmp)
  }
}
df.backbone = as.data.frame(df.backbone)
colnames(df.backbone) = c("data_num","motif_id","input_id","motif_FP_conf","True_conf","False_conf")
rownames(df.backbone) = NULL
# summary results for 4 backbones
df.plot = NULL
for(i in 1:4){
  tmp[[1]] = df.backbone$motif_FP_conf[which(df.backbone$motif_id==2&df.backbone$data_num==i)]
  tmp[[2]] = df.backbone$motif_FP_conf[which(df.backbone$motif_id==3&df.backbone$data_num==i)]
  tmp[[3]] = df.backbone$True_conf[which(df.backbone$motif_id==2&df.backbone$data_num==i)]
  tmp[[4]] = df.backbone$False_conf[which(df.backbone$motif_id==2&df.backbone$data_num==i)]
  for(j in 1:4){
    df.plot = rbind(df.plot,c(i,j,tmp[[j]]))
  }
}
df.plot = as.data.frame(df.plot)
colnames(df.plot) = c("data_num","conf_type","conf")
conf_types = c("Fanout_error","Cascade_error","True","False")
df.plot$conf_type = factor(conf_types[df.plot$conf_type],levels = c("True","Fanout_error","Cascade_error","False"))
backbone_names = c('linear','cycle','bifurcating','converging')
df.plot$backbone = factor(backbone_names[df.plot$data_num],levels = backbone_names)
library(ggplot2)
p = ggplot(df.plot, aes(x = backbone,y = conf, fill = conf_type))  
p = p + labs(title = "Motif errors in 4 backbones") + ylab("Prediction Confidence") + xlab(NULL)
p = p +  geom_bar(stat = "identity", position =position_dodge())
p = p + scale_fill_manual(values=c("#00BFC4","#F8766D","orange","grey"))
p
dev.off()


# effect of parameters X fanout error
dir = 'dyngen/check_controllability_change_all_para/'
data_num = 2

# dataframe of parameters
alphas = 10**seq(-2,1,1)
mhls = c(1,2,5,10)
phls = c(1,2,5,10)
df.para = NULL
for(alpha in alphas){
  for(mhl in mhls){
    for(phl in phls){
      df.para = rbind(df.para,c(alpha,mhl,phl))
    }
  }
}
df.para = as.data.frame(df.para)
colnames(df.para) = c("alpha","mhl","phl")
n_para = nrow(df.para)

alpha = 10
mhl = 1
para_ids = which(df.para$alpha==alpha&df.para$mhl==mhl)
par(mfrow= c(2,2))
df.effect = NULL
for(para_id in para_ids){
  phl = df.para$phl[para_id]
  alpha = df.para$alpha[para_id]
  
  net_id = 1
  rep_id = 0
  seed_id = 1
  if(T){
    load(file = paste(dir,'dataset',data_num,'_net',net_id,"_para",para_id,'.RData',sep = ''))
    
    net = model[["feature_network"]][,c(1,2,5)]
    # modify links: remain only TFs as regulators
    feature_info = model[["feature_info"]]
    TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
    net = net[which(is.element(net$from,TFs)),]
    # remove self regulation
    net = net[which(net$from!=net$to),]
    
    # load NI result 
    if(T){
      
      # load exon result
      load(file = paste(dir,'Network_inference/genie3_output_exon_',data_num,'_net',net_id,"_para",para_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
      exon.mat = weightMatrix
      
      # load inex result
      load(file = paste(dir,'Network_inference/genie3_output_inex_',data_num,'_net',net_id,"_para",para_id,"_",
                        rep_id,"_",seed_id,'.RData',sep = ""))
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
    # mat to linklist
    inex.link = mat_to_link(inex.mat)
    inex.link = inex.link[order(inex.link$weight,decreasing = T),]
    inex.link$confidence = seq(1,0,length.out = nrow(inex.link))
    mat.conf.inex = list_to_mat(inex.link[,c(1,2,4)])
    
    exon.link = mat_to_link(exon.mat)
    exon.link = exon.link[order(exon.link$weight,decreasing = T),]
    exon.link$confidence = seq(1,0,length.out = nrow(exon.link))
    mat.conf.exon = list_to_mat(exon.link[,c(1,2,4)])
    
    # overlap genes
    TFs = rownames(inex.mat)
    targets = colnames(inex.mat)
    net = net[which(is.element(net$from,TFs)),]
    net = net[which(is.element(net$to,targets)),]
    
    
  }
  
  # identify all motifs
  x = graph_from_data_frame(d = net,directed = T)
  z = mat.conf.exon
  motifs = Find_motifs(x)
  print(data_num)
  print(unlist(lapply(motifs,nrow)))
  # print(motifs(x, 3)[c(3,7,5,8)]) 
  
  
  # determine prediction confidence of all instances for motifs
  # par(mfrow= c(2,2))
  motif_names = c("Fan_in","Fan_out","Cascade","FFL")
  input_names = c("inex","exon")
  
  for(motif_id in c(2,3)){
    motif_name = motif_names[motif_id]
    y = motifs[[motif_id]]
    for(input_id in c(1)){
      if(input_id == 1){
        z = mat.conf.inex
      }else{
        z = mat.conf.exon
      }
      input_name = input_names[input_id]
      
      # number of TFs in motif
      num_TFs = apply(y,1,function(y){return(length(intersect(names(V(x)[y]),TFs)))})
      table(num_TFs)
      TF_ids = match(TFs,names(V(x)))
      
      meta_list = NULL
      for(i in 1:3){
        for(j in 1:3){
          if(i!=j){
            meta_list = rbind(meta_list,c(i,j))
          }
        }
      }
      df.conf = NULL
      for(i in 1:nrow(y)){
        tmp = rep(NA,6)
        meta_sub = which(is.element(meta_list[,1],which(is.element(y[i,],TF_ids))))
        for(j in meta_sub){
          tmp[j] = z[names(V(x)[y[i,meta_list[j,1]]]), names(V(x)[y[i,meta_list[j,2]]])]
        }
        df.conf = rbind(df.conf,tmp)
      }
      df.conf = as.data.frame(df.conf)
      colnames(df.conf) = paste(meta_list[,1],meta_list[,2],sep = "")
      rownames(df.conf) = rownames(y)
      
      # determine the background prediction confidence
      true_egdes = net[,1:2]
      true_edges.conf = apply(net,1,function(x){return(z[x[1],x[2]])})
      absent_edges.conf = setdiff(as.numeric(z),true_edges.conf)
      
      # evaluating prediction confidence
      # C12, C13, C23
      if(motif_id==2){
        pvalue = wilcox.test(df.conf$`23`,absent_edges.conf)[["p.value"]]
        title = paste("phl: ",phl," ",input_name," ",motif_name," C23 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`23`,na.rm = T),median(true_edges.conf),median(absent_edges.conf))
      }
      if(motif_id==3){
        pvalue = wilcox.test(df.conf$`13`,absent_edges.conf)[["p.value"]]
        title = paste("phl: ",phl," ",input_name," ",motif_name," C13 pvalue: ",signif(pvalue,2))
        values = c(median(df.conf$`13`),median(true_edges.conf),median(absent_edges.conf))
      }
      boxplot(df.conf$`12`,df.conf$`13`,df.conf$`23`,true_edges.conf,absent_edges.conf,names = c("C12","C13","C23","True","False"),
              main = title)
    }
    
    tmp = c(alpha,mhl,phl,motif_id,input_id,values)
    df.effect = rbind(df.effect,tmp)
  }
}
df.effect = as.data.frame(df.effect)
colnames(df.effect) = c("alpha","mhl","phl","motif_id","input_id","motif_FP_conf","True_conf","False_conf")
rownames(df.effect) = NULL
# plot results
pdf("fanout_error_cycle_phl.pdf",width = 7, height = 7)
if(T){
  par(mfrow  =c(1,1))
  cols = c("red","orange","green","grey")
  id = which(df.effect$motif_id==2)
  plot(df.effect$phl[id],df.effect$motif_FP_conf[id],ylim = c(0,1),xlim = c(0,10),type = "b",col = cols[1],
       main = "Effect of protein half-life on motif errors",xlab = "protein half-life (h)", ylab = "Prediction confidence")
  id = which(df.effect$motif_id==3)
  points(df.effect$phl[id],df.effect$motif_FP_conf[id],type = "b",col = cols[2])
  points(df.effect$phl[id],df.effect$True_conf[id],type = "b",col = cols[3])
  points(df.effect$phl[id],df.effect$False_conf[id],type = "b",col = cols[4])
  legend("bottomleft",c("Fan-out error","Cascade error","True","False"),text.col = cols,col = cols,lty = 1)
}
dev.off()



