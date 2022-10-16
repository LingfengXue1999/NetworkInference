# plot all results     # 2021/9/3
library(ggplot2)

# results for 4 backbones
pdf("./dyngen_results_4backbones.pdf",width = 7,height = 7)
# the variation of configuration for 4 backbones
if(T){
  load("./dyngen/check_controllability/evaluation/variation_of_config_all_backbones.RData")
  df$input = factor(df$input,levels = c("inex","exon"))

  # add random rows
  tmp = df[which(df$input=="exon"),]
  tmp$AUPR = tmp$random_precision
  tmp$input = "random"
  df = rbind(df,tmp)
  df$input = factor(df$input,levels = c("inex","exon","random"))
  # change names
  dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
  df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
  # directly plot
  dataset_list = list('linear','cycle','bifurcating','converging')
  df$data_num = factor(dataset_list[df$data_num],levels = dataset_list)
  p = ggplot(df) + geom_boxplot(aes(x = data_num, y = AUPR,fill = input)) + labs(title = "AUPR for all backbones")
  p = p + ylab("AUPR") + xlab(NULL)
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))
  print(p)
  
  # AUPR ratio plot
  df.ratio = NULL
  for(data_num in 1:4){
    for(net_id in 1:20){
      a = df$AUPR[which(df$data_num==dataset_list[data_num]&df$net_id==net_id&df$input=="premRNA")]
      b = df$AUPR[which(df$data_num==dataset_list[data_num]&df$net_id==net_id&df$input=="mRNA")]
      df.ratio = rbind(df.ratio,c(data_num,net_id,a/b))
    }
  }
  df.ratio = as.data.frame(df.ratio)
  colnames(df.ratio) = c("data_num","net_id","AUPR.ratio")
  # boxplot
  df.ratio$data_num = factor(dataset_list[df.ratio$data_num],levels = dataset_list)
  p = ggplot(df.ratio) + geom_boxplot(aes(x = data_num, y = AUPR.ratio)) + labs(title = "AUPR ratio for all backbones")
  p = p + ylab("AUPR(pre-mRNA) / AUPR(mRNA)") + xlab(NULL)
  p + geom_hline(aes(yintercept=1), colour="black", linetype="dashed")
  print(p)
  
}
dev.off()

# effect of 3 parameters with error bar
pdf(file = "./effect_3para_errorbar.pdf",width = 5,height = 9)
# the variation of configuration: effect of alpha with error bar
if(T){
  load("./dyngen/check_controllability_change_para_error_bar/evaluation/variation_of_config_para_all_backbones.RData")

  alphas = 10**seq(-3,1,4/12)
  alphas = 10**seq(-2,1,1)
  df = df[which(alphas[df$para_id]>=0.01),]
  # add random rows
  tmp = df[which(df$input=="exon"),]
  tmp$AUPR = tmp$random_precision
  tmp$input = "random"
  df = rbind(df,tmp)
  df$input = factor(df$input,levels = c("inex","exon","random"))
  inputs = c("inex","exon","random")
  dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
  inputs = c("premRNA","mRNA","random")
  df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
  
  # calculate mean and sd for 3 reps
  df.sum = NULL
  for(data_num in unique(df$data_num)){
    for(net_id in unique(df$net_id)){
      for(para_id in unique(df$para_id)){
        for(input_id in 1:length(inputs)){
          id = which(df$data_num==data_num&df$net_id==net_id&df$para_id==para_id&df$input==inputs[input_id])
          tmp = c(data_num,net_id,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
          df.sum = rbind(df.sum,tmp)
        }
      }
    }
  }
  
  df.sum = as.data.frame(df.sum)
  colnames(df.sum) = c("data_num","net_id","para_id","input_id","AUPR.m","AUPR.sd")
  df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
  
  
  # further summary: summary all configurations
  df.sum2 = NULL
  for(data_num in unique(df$data_num)){
    for(para_id in unique(df.sum$para_id)){
      for(input_id in 1:length(inputs)){
        id = which(df.sum$data_num==data_num&df.sum$para_id==para_id&df.sum$input==inputs[input_id])
        tmp = c(data_num,para_id,input_id,mean(df.sum$AUPR.m[id]),sd(df.sum$AUPR.m[id]))
        df.sum2 = rbind(df.sum2,tmp)
      }
    }
  }
  df.sum2 = as.data.frame(df.sum2)
  colnames(df.sum2) = c("data_num","para_id","input_id","AUPR.m.m","AUPR.m.sd")
  df.sum2$input = factor(inputs[df.sum2$input_id],levels = inputs)
  # error bar plot for all configurations
  l.plot = NULL
  for(i in 1:length(unique(df$data_num))){
    df.plot = df.sum2[which(df.sum2$data_num==i),]
    df.plot$alpha = (alphas[df.plot$para_id])
    p = ggplot(df.plot,aes(x = alpha, y = AUPR.m.m,colour = input, group = input, fill = input)) + 
      geom_errorbar(aes(ymin = AUPR.m.m - AUPR.m.sd, ymax = AUPR.m.m + AUPR.m.sd), width = 0.2, position = position_dodge(0.1))
    p = p + geom_line(position = position_dodge(0.1), size = 1)
    p = p + scale_x_log10()+ labs(title = paste("Backbone",i)) + ylab("AUPR")
    p = p + scale_color_manual(values=c("#F8766D","#00BFC4","grey"))
    
    l.plot[[i]] = p
  }
  library(cowplot)
  i = 0
  print(  plot_grid(l.plot[[i*4+1]],l.plot[[i*4+2]],l.plot[[i*4+3]],l.plot[[i*4+4]],ncol = 1,align = "v"))
  
  # summary results for 4 backbones: AUPR (max value) / AUPR (min value)
  df.plot = df.sum2[which(df.sum2$para_id==1),]
  df.plot$ratio = df.sum2$AUPR.m.m[which(df.sum2$para_id==4)]/df.sum2$AUPR.m.m[which(df.sum2$para_id==1)]
  df.plot = df.plot[which(df.plot$input_id==1|df.plot$input_id==2),]
  df.plot$input = factor(df.plot$input, levels = c("premRNA","mRNA"))
  # plot AUPR
  p = ggplot(df.plot, aes(x = data_num,y = ratio, fill = input))  + labs(title = "Effect of transcription rate")
  p = p + ylab("AUPR ratio") + xlab(NULL)
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4")) 
  p = p +  geom_bar(stat = "identity", position =position_dodge())
  p = p + geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")
  p.1 = p
}

# effect of mRNA half-life
if(T){
  load("./dyngen/check_controllability_change_para_error_bar_mhl/evaluation/variation_of_config_para_all_backbones.RData")
  
  mhls = c(1,2,5,10)
  # add random rows
  tmp = df[which(df$input=="exon"),]
  tmp$AUPR = tmp$random_precision
  tmp$input = "random"
  df = rbind(df,tmp)
  df$input = factor(df$input,levels = c("inex","exon","random"))
  inputs = c("inex","exon","random")
  dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
  inputs = c("premRNA","mRNA","random")
  df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
  
  # calculate mean and sd for 3 reps
  df.sum = NULL
  for(data_num in unique(df$data_num)){
    for(net_id in unique(df$net_id)){
      for(para_id in unique(df$para_id)){
        for(input_id in 1:length(inputs)){
          id = which(df$data_num==data_num&df$net_id==net_id&df$para_id==para_id&df$input==inputs[input_id])
          tmp = c(data_num,net_id,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
          df.sum = rbind(df.sum,tmp)
        }
      }
    }
  }
  
  df.sum = as.data.frame(df.sum)
  colnames(df.sum) = c("data_num","net_id","para_id","input_id","AUPR.m","AUPR.sd")
  df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
  
  
  # further summary: summary all configurations
  df.sum2 = NULL
  for(data_num in unique(df$data_num)){
    for(para_id in unique(df.sum$para_id)){
      for(input_id in 1:length(inputs)){
        id = which(df.sum$data_num==data_num&df.sum$para_id==para_id&df.sum$input==inputs[input_id])
        tmp = c(data_num,para_id,input_id,mean(df.sum$AUPR.m[id]),sd(df.sum$AUPR.m[id]))
        df.sum2 = rbind(df.sum2,tmp)
      }
    }
  }
  df.sum2 = as.data.frame(df.sum2)
  colnames(df.sum2) = c("data_num","para_id","input_id","AUPR.m.m","AUPR.m.sd")
  df.sum2$input = factor(inputs[df.sum2$input_id],levels = inputs)
  # error bar plot for all configurations
  l.plot = NULL
  for(i in 1:length(unique(df$data_num))){
    df.plot = df.sum2[which(df.sum2$data_num==i),]
    df.plot$mhl = (mhls[df.plot$para_id])
    p = ggplot(df.plot,aes(x = mhl, y = AUPR.m.m,colour = input, group = input, fill = input)) + 
      geom_errorbar(aes(ymin = AUPR.m.m - AUPR.m.sd, ymax = AUPR.m.m + AUPR.m.sd), width = 0.2, position = position_dodge(0.1))
    p = p + geom_line(position = position_dodge(0.1), size = 1)
    p = p + labs(title = paste("Backbone",i)) + ylab("AUPR") + xlab("mRNA half-life (h)")
    p = p + scale_color_manual(values=c("#F8766D","#00BFC4","grey"))
    
    l.plot[[i]] = p
  }
  library(cowplot)
  i = 0
  print(  plot_grid(l.plot[[i*4+1]],l.plot[[i*4+2]],l.plot[[i*4+3]],l.plot[[i*4+4]],ncol = 1,align = "v"))
  
  # summary results for 4 backbones: AUPR (max value) / AUPR (min value)
  df.plot = df.sum2[which(df.sum2$para_id==1),]
  df.plot$ratio = df.sum2$AUPR.m.m[which(df.sum2$para_id==4)]/df.sum2$AUPR.m.m[which(df.sum2$para_id==1)]
  df.plot = df.plot[which(df.plot$input_id==1|df.plot$input_id==2),]
  df.plot$input = factor(df.plot$input, levels = c("premRNA","mRNA"))
  # plot AUPR
  p = ggplot(df.plot, aes(x = data_num,y = ratio, fill = input))  + labs(title = "Effect of mRNA half-life")
  p = p + ylab("AUPR ratio") + xlab(NULL)
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4")) 
  p = p +  geom_bar(stat = "identity", position =position_dodge())
  p = p + geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")
  p.2 = p
  
  
  
}

# effect of protein half-life of TF
if(T){
  load("./dyngen/check_controllability_change_para_error_bar_phl_TF/evaluation/variation_of_config_para_all_backbones.RData")
  
  phls = c(1,2,5,10)
  # add random rows
  tmp = df[which(df$input=="exon"),]
  tmp$AUPR = tmp$random_precision
  tmp$input = "random"
  df = rbind(df,tmp)
  df$input = factor(df$input,levels = c("inex","exon","random"))
  inputs = c("inex","exon","random")
  dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
  inputs = c("premRNA","mRNA","random")
  df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
  
  # calculate mean and sd for 3 reps
  df.sum = NULL
  for(data_num in unique(df$data_num)){
    for(net_id in unique(df$net_id)){
      for(para_id in unique(df$para_id)){
        for(input_id in 1:length(inputs)){
          id = which(df$data_num==data_num&df$net_id==net_id&df$para_id==para_id&df$input==inputs[input_id])
          tmp = c(data_num,net_id,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
          df.sum = rbind(df.sum,tmp)
        }
      }
    }
  }
  
  df.sum = as.data.frame(df.sum)
  colnames(df.sum) = c("data_num","net_id","para_id","input_id","AUPR.m","AUPR.sd")
  df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
  
  
  # further summary: summary all configurations
  df.sum2 = NULL
  for(data_num in unique(df$data_num)){
    for(para_id in unique(df.sum$para_id)){
      for(input_id in 1:length(inputs)){
        id = which(df.sum$data_num==data_num&df.sum$para_id==para_id&df.sum$input==inputs[input_id])
        tmp = c(data_num,para_id,input_id,mean(df.sum$AUPR.m[id]),sd(df.sum$AUPR.m[id]))
        df.sum2 = rbind(df.sum2,tmp)
      }
    }
  }
  df.sum2 = as.data.frame(df.sum2)
  colnames(df.sum2) = c("data_num","para_id","input_id","AUPR.m.m","AUPR.m.sd")
  df.sum2$input = factor(inputs[df.sum2$input_id],levels = inputs)
  # error bar plot for all configurations
  l.plot = NULL
  for(i in 1:length(unique(df$data_num))){
    df.plot = df.sum2[which(df.sum2$data_num==i),]
    df.plot$phl = (phls[df.plot$para_id])
    p = ggplot(df.plot,aes(x = phl, y = AUPR.m.m,colour = input, group = input, fill = input)) + 
      geom_errorbar(aes(ymin = AUPR.m.m - AUPR.m.sd, ymax = AUPR.m.m + AUPR.m.sd), width = 0.2, position = position_dodge(0.1))
    p = p + geom_line(position = position_dodge(0.1), size = 1)
    p = p + labs(title = paste("Backbone",i)) + ylab("AUPR") + xlab("protein half-life of TFs (h)")
    p = p + scale_color_manual(values=c("#F8766D","#00BFC4","grey"))
    
    l.plot[[i]] = p
  }
  library(cowplot)
  i = 0
  print(  plot_grid(l.plot[[i*4+1]],l.plot[[i*4+2]],l.plot[[i*4+3]],l.plot[[i*4+4]],ncol = 1,align = "v"))
  
  # summary results for 4 backbones: AUPR (max value) / AUPR (min value)
  df.plot = df.sum2[which(df.sum2$para_id==1),]
  df.plot$ratio = df.sum2$AUPR.m.m[which(df.sum2$para_id==4)]/df.sum2$AUPR.m.m[which(df.sum2$para_id==1)]
  df.plot = df.plot[which(df.plot$input_id==1|df.plot$input_id==2),]
  df.plot$input = factor(df.plot$input, levels = c("premRNA","mRNA"))
  # plot AUPR
  p = ggplot(df.plot, aes(x = data_num,y = ratio, fill = input))  + labs(title = "Effect of protein half-life")
  p = p + ylab("AUPR ratio") + xlab("backbone")
  p = p + scale_fill_manual(values=c("#F8766D","#00BFC4")) 
  p = p +  geom_bar(stat = "identity", position =position_dodge())
  p = p + geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")
  p.3 = p
  
  
  
}
dev.off()

# summary effect of 3 parameters 
pdf(file = "./summary_effect_3para.pdf",width = 9,height = 9)
if(T){
  library(cowplot)
  plot_grid(p.1,p.2,p.3,ncol = 1,align = "v")
}
dev.off()

# effect of 3 parameters with error bar -- only one backbone
pdf(file = "./effect_3para_errorbar_converging.pdf",width = 6,height = 4)
if(T){
  # the variation of configuration: effect of alpha with error bar
  if(T){
    load("./dyngen/check_controllability_change_para_error_bar/evaluation/variation_of_config_para_all_backbones.RData")
    
    alphas = 10**seq(-3,1,4/12)
    alphas = 10**seq(-2,1,1)
    df = df[which(alphas[df$para_id]>=0.01),]
    # add random rows
    tmp = df[which(df$input=="exon"),]
    tmp$AUPR = tmp$random_precision
    tmp$input = "random"
    df = rbind(df,tmp)
    df$input = factor(df$input,levels = c("inex","exon","random"))
    inputs = c("inex","exon","random")
    dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
    inputs = c("premRNA","mRNA","random")
    df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
    
    # calculate mean and sd for 3 reps
    df.sum = NULL
    for(data_num in unique(df$data_num)){
      for(net_id in unique(df$net_id)){
        for(para_id in unique(df$para_id)){
          for(input_id in 1:length(inputs)){
            id = which(df$data_num==data_num&df$net_id==net_id&df$para_id==para_id&df$input==inputs[input_id])
            tmp = c(data_num,net_id,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
            df.sum = rbind(df.sum,tmp)
          }
        }
      }
    }
    
    df.sum = as.data.frame(df.sum)
    colnames(df.sum) = c("data_num","net_id","para_id","input_id","AUPR.m","AUPR.sd")
    df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
    
    
    # further summary: summary all configurations
    df.sum2 = NULL
    for(data_num in unique(df$data_num)){
      for(para_id in unique(df.sum$para_id)){
        for(input_id in 1:length(inputs)){
          id = which(df.sum$data_num==data_num&df.sum$para_id==para_id&df.sum$input==inputs[input_id])
          tmp = c(data_num,para_id,input_id,mean(df.sum$AUPR.m[id]),sd(df.sum$AUPR.m[id]))
          df.sum2 = rbind(df.sum2,tmp)
        }
      }
    }
    df.sum2 = as.data.frame(df.sum2)
    colnames(df.sum2) = c("data_num","para_id","input_id","AUPR.m.m","AUPR.m.sd")
    df.sum2$input = factor(inputs[df.sum2$input_id],levels = inputs)
    # error bar plot for all configurations
    l.plot = NULL
    for(i in 1:length(unique(df$data_num))){
      df.plot = df.sum2[which(df.sum2$data_num==i),]
      df.plot$alpha = (alphas[df.plot$para_id])
      p = ggplot(df.plot,aes(x = alpha, y = AUPR.m.m,colour = input, group = input, fill = input)) + 
        geom_errorbar(aes(ymin = AUPR.m.m - AUPR.m.sd, ymax = AUPR.m.m + AUPR.m.sd), width = 0.2, position = position_dodge(0.1))
      p = p + geom_line(position = position_dodge(0.1), size = 1)
      p = p + scale_x_log10()+ labs(title = paste("Backbone",i)) + ylab("AUPR")
      p = p + scale_color_manual(values=c("#F8766D","#00BFC4","grey"))
      
      l.plot[[i]] = p
    }
    library(cowplot)
    i = 0
    print(  l.plot[[i*4+4]])
    
    # summary results for 4 backbones: AUPR (max value) / AUPR (min value)
    df.plot = df.sum2[which(df.sum2$para_id==1),]
    df.plot$ratio = df.sum2$AUPR.m.m[which(df.sum2$para_id==4)]/df.sum2$AUPR.m.m[which(df.sum2$para_id==1)]
    df.plot = df.plot[which(df.plot$input_id==1|df.plot$input_id==2),]
    df.plot$input = factor(df.plot$input, levels = c("premRNA","mRNA"))
    # plot AUPR
    p = ggplot(df.plot, aes(x = data_num,y = ratio, fill = input))  + labs(title = "Effect of transcription rate")
    p = p + ylab("AUPR ratio") + xlab(NULL)
    p = p + scale_fill_manual(values=c("#F8766D","#00BFC4")) 
    p = p +  geom_bar(stat = "identity", position =position_dodge())
    p = p + geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")
    p.1 = p
  }
  
  # effect of mRNA half-life
  if(T){
    load("./dyngen/check_controllability_change_para_error_bar_mhl/evaluation/variation_of_config_para_all_backbones.RData")
    
    mhls = c(1,2,5,10)
    # add random rows
    tmp = df[which(df$input=="exon"),]
    tmp$AUPR = tmp$random_precision
    tmp$input = "random"
    df = rbind(df,tmp)
    df$input = factor(df$input,levels = c("inex","exon","random"))
    inputs = c("inex","exon","random")
    dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
    inputs = c("premRNA","mRNA","random")
    df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
    
    # calculate mean and sd for 3 reps
    df.sum = NULL
    for(data_num in unique(df$data_num)){
      for(net_id in unique(df$net_id)){
        for(para_id in unique(df$para_id)){
          for(input_id in 1:length(inputs)){
            id = which(df$data_num==data_num&df$net_id==net_id&df$para_id==para_id&df$input==inputs[input_id])
            tmp = c(data_num,net_id,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
            df.sum = rbind(df.sum,tmp)
          }
        }
      }
    }
    
    df.sum = as.data.frame(df.sum)
    colnames(df.sum) = c("data_num","net_id","para_id","input_id","AUPR.m","AUPR.sd")
    df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
    
    
    # further summary: summary all configurations
    df.sum2 = NULL
    for(data_num in unique(df$data_num)){
      for(para_id in unique(df.sum$para_id)){
        for(input_id in 1:length(inputs)){
          id = which(df.sum$data_num==data_num&df.sum$para_id==para_id&df.sum$input==inputs[input_id])
          tmp = c(data_num,para_id,input_id,mean(df.sum$AUPR.m[id]),sd(df.sum$AUPR.m[id]))
          df.sum2 = rbind(df.sum2,tmp)
        }
      }
    }
    df.sum2 = as.data.frame(df.sum2)
    colnames(df.sum2) = c("data_num","para_id","input_id","AUPR.m.m","AUPR.m.sd")
    df.sum2$input = factor(inputs[df.sum2$input_id],levels = inputs)
    # error bar plot for all configurations
    l.plot = NULL
    for(i in 1:length(unique(df$data_num))){
      df.plot = df.sum2[which(df.sum2$data_num==i),]
      df.plot$mhl = (mhls[df.plot$para_id])
      p = ggplot(df.plot,aes(x = mhl, y = AUPR.m.m,colour = input, group = input, fill = input)) + 
        geom_errorbar(aes(ymin = AUPR.m.m - AUPR.m.sd, ymax = AUPR.m.m + AUPR.m.sd), width = 0.2, position = position_dodge(0.1))
      p = p + geom_line(position = position_dodge(0.1), size = 1)
      p = p + labs(title = paste("Backbone",i)) + ylab("AUPR") + xlab("mRNA half-life (h)")
      p = p + scale_color_manual(values=c("#F8766D","#00BFC4","grey"))
      
      l.plot[[i]] = p
    }
    library(cowplot)
    i = 0
    print(l.plot[[i*4+4]])
    
    # summary results for 4 backbones: AUPR (max value) / AUPR (min value)
    df.plot = df.sum2[which(df.sum2$para_id==1),]
    df.plot$ratio = df.sum2$AUPR.m.m[which(df.sum2$para_id==4)]/df.sum2$AUPR.m.m[which(df.sum2$para_id==1)]
    df.plot = df.plot[which(df.plot$input_id==1|df.plot$input_id==2),]
    df.plot$input = factor(df.plot$input, levels = c("premRNA","mRNA"))
    # plot AUPR
    p = ggplot(df.plot, aes(x = data_num,y = ratio, fill = input))  + labs(title = "Effect of mRNA half-life")
    p = p + ylab("AUPR ratio") + xlab(NULL)
    p = p + scale_fill_manual(values=c("#F8766D","#00BFC4")) 
    p = p +  geom_bar(stat = "identity", position =position_dodge())
    p = p + geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")
    p.2 = p
    
    
    
  }
  
  # effect of protein half-life of TF
  if(T){
    load("./dyngen/check_controllability_change_para_error_bar_phl_TF/evaluation/variation_of_config_para_all_backbones.RData")
    
    phls = c(1,2,5,10)
    # add random rows
    tmp = df[which(df$input=="exon"),]
    tmp$AUPR = tmp$random_precision
    tmp$input = "random"
    df = rbind(df,tmp)
    df$input = factor(df$input,levels = c("inex","exon","random"))
    inputs = c("inex","exon","random")
    dict = data.frame(inex = "premRNA",exon = "mRNA", random = "random")
    inputs = c("premRNA","mRNA","random")
    df$input = factor(as.character(as.matrix(dict[df$input])),levels = as.character(as.matrix(dict)) )
    
    # calculate mean and sd for 3 reps
    df.sum = NULL
    for(data_num in unique(df$data_num)){
      for(net_id in unique(df$net_id)){
        for(para_id in unique(df$para_id)){
          for(input_id in 1:length(inputs)){
            id = which(df$data_num==data_num&df$net_id==net_id&df$para_id==para_id&df$input==inputs[input_id])
            tmp = c(data_num,net_id,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
            df.sum = rbind(df.sum,tmp)
          }
        }
      }
    }
    
    df.sum = as.data.frame(df.sum)
    colnames(df.sum) = c("data_num","net_id","para_id","input_id","AUPR.m","AUPR.sd")
    df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
    
    
    # further summary: summary all configurations
    df.sum2 = NULL
    for(data_num in unique(df$data_num)){
      for(para_id in unique(df.sum$para_id)){
        for(input_id in 1:length(inputs)){
          id = which(df.sum$data_num==data_num&df.sum$para_id==para_id&df.sum$input==inputs[input_id])
          tmp = c(data_num,para_id,input_id,mean(df.sum$AUPR.m[id]),sd(df.sum$AUPR.m[id]))
          df.sum2 = rbind(df.sum2,tmp)
        }
      }
    }
    df.sum2 = as.data.frame(df.sum2)
    colnames(df.sum2) = c("data_num","para_id","input_id","AUPR.m.m","AUPR.m.sd")
    df.sum2$input = factor(inputs[df.sum2$input_id],levels = inputs)
    # error bar plot for all configurations
    l.plot = NULL
    for(i in 1:length(unique(df$data_num))){
      df.plot = df.sum2[which(df.sum2$data_num==i),]
      df.plot$phl = (phls[df.plot$para_id])
      p = ggplot(df.plot,aes(x = phl, y = AUPR.m.m,colour = input, group = input, fill = input)) + 
        geom_errorbar(aes(ymin = AUPR.m.m - AUPR.m.sd, ymax = AUPR.m.m + AUPR.m.sd), width = 0.2, position = position_dodge(0.1))
      p = p + geom_line(position = position_dodge(0.1), size = 1)
      p = p + labs(title = paste("Backbone",i)) + ylab("AUPR") + xlab("protein half-life of TFs (h)")
      p = p + scale_color_manual(values=c("#F8766D","#00BFC4","grey"))
      
      l.plot[[i]] = p
    }
    library(cowplot)
    i = 0
    print(  l.plot[[i*4+4]])
    
    # summary results for 4 backbones: AUPR (max value) / AUPR (min value)
    df.plot = df.sum2[which(df.sum2$para_id==1),]
    df.plot$ratio = df.sum2$AUPR.m.m[which(df.sum2$para_id==4)]/df.sum2$AUPR.m.m[which(df.sum2$para_id==1)]
    df.plot = df.plot[which(df.plot$input_id==1|df.plot$input_id==2),]
    df.plot$input = factor(df.plot$input, levels = c("premRNA","mRNA"))
    # plot AUPR
    p = ggplot(df.plot, aes(x = data_num,y = ratio, fill = input))  + labs(title = "Effect of protein half-life")
    p = p + ylab("AUPR ratio") + xlab("backbone")
    p = p + scale_fill_manual(values=c("#F8766D","#00BFC4")) 
    p = p +  geom_bar(stat = "identity", position =position_dodge())
    p = p + geom_hline(aes(yintercept = 1),colour = "black",linetype = "dashed")
    p.3 = p
    
    
    
  }
}
dev.off()




# synergic effect of parameters
pdf(file = "./synergic_effect_3para.pdf",width = 5,height = 5)
data_num=1
for(data_num in 1:2){
  if(data_num == 1){
    load("./dyngen/check_controllability_change_all_para/evaluation/variation_of_config_para_all_backbones.RData")
  }else if(data_num ==2){
    load("./dyngen/check_controllability_change_all_para/evaluation/variation_of_config_para_all_backbones_2.RData")
  }
  
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
  
  # add random rows
  tmp = df[which(df$input=="exon"),]
  tmp$AUPR = tmp$random_precision
  tmp$input = "random"
  df = rbind(df,tmp)
  df$input = factor(df$input,levels = c("inex","exon","random"))
  inputs = c("inex","exon","random")
  
  # combine metadata
  df = cbind(df,df.para[df$para_id,])

  # calculate mean and sd for configs
  df.sum = NULL
  for(data_num in unique(df$data_num)){
    for(para_id in unique(df$para_id)){
      for(input_id in 1:length(inputs)){
        id = which(df$data_num==data_num&df$para_id==para_id&df$input==inputs[input_id])
        tmp = c(data_num,para_id,input_id,mean(df$AUPR[id]),sd(df$AUPR[id]))
        df.sum = rbind(df.sum,tmp)
      }
    }
  }
  df.sum = as.data.frame(df.sum)
  colnames(df.sum) = c("data_num","para_id","input_id","AUPR.m","AUPR.sd")
  df.sum$input = factor(inputs[df.sum$input_id],levels = inputs)
  
  # summary data
  inex.m = df.sum[which(df.sum$input_id==1),]
  exon.m = df.sum[which(df.sum$input_id==2),]
  rand.m = df.sum[which(df.sum$input_id==3),]
  df.m = data.frame(inex.m = inex.m$AUPR.m, exon.m = exon.m$AUPR.m, rand.m = rand.m$AUPR.m)
  df.m = cbind(df.m,df.para)
  
  backbone_names = c("linear","cycle","bifurcating","converging")
  # 2d heatmap plot
  if(T){
    # two half-life plot -- exon results
    df.plot = df.m[which(df.m$alpha==10),]
    p = ggplot(df.plot, aes(x=mhl, y=phl,color = exon.m)) + geom_point(size = 10)
    p = p +labs(title=paste(data_num,"Effect of two parameters on mRNA"),colour = "AUPR")+ xlab("mRNA half-life (h)") +ylab("protein half-life (h)")
    p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
    print(p)
    
    # alpha - protien half-life plot
    df.plot = df.m[which(df.m$mhl==10),]
    df.plot$diff = df.plot$inex.m - df.plot$exon.m
    p = ggplot(df.plot, aes(x=alpha, y=phl,color = diff)) + geom_point(size = 10) 
    p = p + xlab("alpha (1/min)") + ylab("protein half-life (h)") + labs(title=paste(data_num,"Compare pre-mRNA with mRNA"))
    p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red") + scale_x_log10()
    print(p)
  }
  
  if(data_num ==2){
    # effect of alpha; mhl, phl as parameters
    par(mfrow = c(1,1))
    mhl = 10
    phl = 1
    id = which(df.m$mhl==mhl&df.m$phl==phl)
    df.sub = df.m[id,]
    cols1 = rgb(255, 0, 0, 2**seq(log2(255),log2(255/5),length.out = 4), maxColorValue=255)
    
    
    plot(df.sub$alpha,df.sub$inex.m,col = "red",main = "Cycle: protein half-life X alpha",ylim = c(0,1),ylab = "AUPR",xlab = "alpha",
         log = "x",type = "b",lty = 2)
    
    for(phl_id in 1:length(unique(df.m$phl))){
      phl = unique(df.m$phl)[phl_id]
      id = which(df.m$mhl==mhl&df.m$phl==phl)
      df.sub = df.m[id,]
      lines(df.sub$alpha,df.sub$inex.m,col = cols1[phl_id],pch = phl_id,type = "b",lty = 2)
    }
    legend("bottomleft",paste("phl=",c(1,2,5,10),"h",sep = ""),text.col = cols1,col = cols1,lty = 2,pch = c(1,2,3,4))
  }
  
  
}
dev.off()
