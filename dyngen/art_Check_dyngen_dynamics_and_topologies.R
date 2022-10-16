# Check dyngen dynamics  # 2021/9/2
# function to detect dynamics
# remove the effect of fluctuation (var(X) = mean(X))
det_dyn = function(x){
  thr = (max(x) + min(x))/2
  sd = sqrt(thr)
  # check >thr or <thr for each point
  i.thr = NULL
  state = x[1]>thr
  for(i in 1:length(x)){
    state.new = x[i]>thr
    if(state!=state.new){
      i.thr = c(i.thr,i)
    }
    state = state.new
  }
  if(length(i.thr)>=2){
    return("oscillation")
  }else if(length(i.thr)>0){
    if(x[1]>thr){
      return("down")
    }else{
      return("up")
    }
  }else{
    return("constant")
  }
}
# dynamics for 4 backbones
dataset_names = list('linear','cycle','bifurcating','converging')
data_num = 4
TF_names_list = list()
for(data_num in 1:4){
  dir = 'dyngen/check_controllability/'
  net_id = 1
  
  # load input matrix
  load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
  simu_counts = model[["simulations"]][["counts"]]
  t = model[["simulations"]][["meta"]][["sim_time"]]
  
  # count data
  n_genes = ncol(simu_counts)/3
  genes = model[["feature_info"]]$feature_id
  s = t(as.matrix(simu_counts[,(1:n_genes)+n_genes]))
  u = t(as.matrix(simu_counts[,1:n_genes]))
  rownames(s) = rownames(u) = genes
  
  # time data
  t.per = which(t==t[1])[2]-1    # extract period of simulation time
  n_t = length(t)/t.per
  
  net = model[["feature_network"]]
  tmp = intersect(1:t.per,which(t>=0))
  
  # investigate the dynamics of each TF
  feature_info = model[["feature_info"]]
  TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
  par(mfrow = c(2,2))
  for(i in 1:length(TFs)){
    TF = TFs[i]
    x = s[TF,tmp]
    plot(t[tmp],x,type = "l",ylim = c(0,max(x)),main = TF)
    thr = (min(x)+max(x))/2
    abline(h = thr)
  }
  
  TF_names_list[[length(TF_names_list)+1]]= TFs
}
# 4 types of dynamics
pdf("./4_types_of_dynamics.pdf",width = 9,height = 9)
par(mfrow = c(2,2))
if(T){
  data_num = 1
  dir = 'dyngen/check_controllability/'
  net_id = 1
  
  # load input matrix
  load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
  simu_counts = model[["simulations"]][["counts"]]
  t = model[["simulations"]][["meta"]][["sim_time"]]
  
  # count data
  n_genes = ncol(simu_counts)/3
  genes = model[["feature_info"]]$feature_id
  s = t(as.matrix(simu_counts[,(1:n_genes)+n_genes]))
  u = t(as.matrix(simu_counts[,1:n_genes]))
  rownames(s) = rownames(u) = genes
  
  # time data
  t.per = which(t==t[1])[2]-1    # extract period of simulation time
  n_t = length(t)/t.per
  
  net = model[["feature_network"]]
  tmp = intersect(1:t.per,which(t>=0))
  
  plot(t[tmp],s["M2_TF1",tmp],type = "l",ylim = c(0,max(s["M2_TF1",tmp])),main = "Up dynamics",xlab = "t (h)",ylab = "mRNA level")
  plot(t[tmp],s["M3_TF1",tmp],type = "l",ylim = c(0,max(s["M3_TF1",tmp])),main = "Down dynamics",xlab = "t (h)",ylab = "mRNA level")
  
  data_num = 2
  dir = 'dyngen/check_controllability/'
  net_id = 1
  
  # load input matrix
  load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
  simu_counts = model[["simulations"]][["counts"]]
  t = model[["simulations"]][["meta"]][["sim_time"]]
  
  # count data
  n_genes = ncol(simu_counts)/3
  genes = model[["feature_info"]]$feature_id
  s = t(as.matrix(simu_counts[,(1:n_genes)+n_genes]))
  u = t(as.matrix(simu_counts[,1:n_genes]))
  rownames(s) = rownames(u) = genes
  
  # time data
  t.per = which(t==t[1])[2]-1    # extract period of simulation time
  n_t = length(t)/t.per
  
  net = model[["feature_network"]]
  tmp = intersect(1:t.per,which(t>=0))
  
  plot(t[tmp],s["M1_TF1",tmp],type = "l",ylim = c(0,max(s["M1_TF1",tmp])),main = "Oscillation dynamics",xlab = "t (h)",ylab = "mRNA level")
  
  data_num = 1
  dir = 'dyngen/check_controllability/'
  net_id = 1
  
  # load input matrix
  load(file = paste(dir,'dataset',data_num,'_net',net_id,'.RData',sep = ''))
  simu_counts = model[["simulations"]][["counts"]]
  t = model[["simulations"]][["meta"]][["sim_time"]]
  
  # count data
  n_genes = ncol(simu_counts)/3
  genes = model[["feature_info"]]$feature_id
  s = t(as.matrix(simu_counts[,(1:n_genes)+n_genes]))
  u = t(as.matrix(simu_counts[,1:n_genes]))
  rownames(s) = rownames(u) = genes
  
  # time data
  t.per = which(t==t[1])[2]-1    # extract period of simulation time
  n_t = length(t)/t.per
  
  net = model[["feature_network"]]
  tmp = intersect(1:t.per,which(t>=0))
  
  plot(t[tmp],s["M1_TF1",tmp],type = "l",ylim = c(0,max(s["M1_TF1",tmp])),main = "Constant dynamics",xlab = "t (h)",ylab = "mRNA level")
  
}
dev.off()

# results for dynamics of 4 backbones
pdf("./dynamics_of_4backbones.pdf",width = 16,height = 9)
if(T){
  dynamics_type = c("up","down","oscillation","constant")
  TF_dynamics_list = list()
  TF_dynamics_list[[1]] = dynamics_type[c(4,1,2,1,1)]
  TF_dynamics_list[[2]] = dynamics_type[c(3,3,3,3,3)]
  TF_dynamics_list[[3]] = dynamics_type[c(4,4,4,4,4,1,1,1,1,1,1,3,2,3,1,4,1,1,1,1,4,4,4,1,1,2,1,4,4,4,4)]
  TF_dynamics_list[[4]] = dynamics_type[c(2,4,4,2,3,1,4,4,1,1,2,1,1,1,1,1)]
  
  # pie plot
  i = 1
  
  par(mfrow = c(1,4))
  for(i in 1:4){
    TF_dynamics_list[[i]] = factor(TF_dynamics_list[[i]],levels = c("up","down","oscillation","constant"))
    info = table(TF_dynamics_list[[i]])
    cols= c("red","blue","purple","grey")
    pie(info, labels = info, col=cols,main = paste(dataset_names[i]))
    legend("topright",c("up","down","oscillation","constant"),text.col = cols,col = cols)
    
  }
}
dev.off()

# dynamics for fan-out error in cycle backbone -- control 3 parameters
pdf("./cycle_fanout_error_dynamics.pdf",width = 9,height = 9)
if(T){
  dir = 'dyngen/check_controllability_change_all_para/'
  data_num = 2
  net_id = 1
  para_id = 1
  # parameters
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
  
  phl = 10
  mhl = 10
  para_ids = which(df.para$mhl==mhl&df.para$phl==phl)
  
  # View dynamics
  i = 2
  para_id = para_ids[i]
  # load input matrix
  load(file = paste(dir,'dataset',data_num,'_net',net_id,"_para",para_id,'.RData',sep = ''))
  simu_counts = model[["simulations"]][["counts"]]
  t = model[["simulations"]][["meta"]][["sim_time"]]
  
  # count data
  n_genes = ncol(simu_counts)/3
  genes = model[["feature_info"]]$feature_id
  s = t(as.matrix(simu_counts[,(1:n_genes)+n_genes]))
  u = t(as.matrix(simu_counts[,1:n_genes]))
  rownames(s) = rownames(u) = genes
  
  # time data
  t.per = which(t==t[1])[2]-1    # extract period of simulation time
  n_t = length(t)/t.per
  
  net = model[["feature_network"]]
  tmp = intersect(1:t.per,which(t>=0))
  
  TF1 = "M4_TF1"
  TF2 = "M5_TF1"
  Target ="Target18"
  length(which(net$from==TF1&net$to==TF2))
  length(which(net$from==TF1&net$to==Target))
  length(which(net$from==TF2&net$to==Target))
  
  
  plot(t[tmp],s[TF1,tmp],type = "l",col = "red",ylim = c(0,200),xlab = "t (h)", 
       ylab = "mRNA level",main = "TF1 -> TF2, TF1 -> Target")
  points(t[tmp],s[TF2,tmp],type = "l",col = "blue")
  points(t[tmp],s[Target,tmp],type = "l",col = "black")
  cols = c("red","blue","black")
  legend("topright",c("TF1","TF2","Target"),text.col = cols,col = cols,lty = 1)
}
dev.off()

# Check network topology of dyngen # 2021/8/21
library(tidyverse)
library(dyngen)
library("ggplot2")
library(igraph)
library(Matrix)

dataset_list = list('linear_simple','cycle_simple','bifurcating','converging')
backbone_names = c("linear","cycle","bifurcating","converging")

# plot 4 backbones
pdf("./backbone_topology_of_dyngen.pdf",width = 16,height = 9)
par(mfrow = c(1,4))
# determin backbone
for(data_num in 1:4){
  set.seed(0) # set backbone for bifurcating 
  backbone = list_backbones()[[dataset_list[[data_num]]]]()
  
  net = backbone$module_network
  x = graph_from_data_frame(d = net,directed = T)
  E(x)$effect = factor(E(x)$effect,levels = c(1,-1))
  cols = c("red","blue")
  E(x)$color = cols[E(x)$effect]
  E(x)$arrow.size = 0.3
  V(x)$color = "white"
  V(x)$size = 8
  if(data_num==1|data_num==2){
    l = layout_as_star(x)
  }
  if(data_num==3|data_num==4){
    l = layout.reingold.tilford(x)
  }
  
  # l = layout_with_drl(x)
  plot(x,layout = l,main= backbone_names[data_num],vertex.label=NA)
  cols = c("red","blue")
  legend("topleft",c("Activation","Inhibition"),text.col = cols,col = cols,lty = 1)
}
dev.off()
