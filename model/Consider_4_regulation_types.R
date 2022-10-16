# Consider 4 types of regulation   # 2021/4/30
# analytical solution for unspliced mRNA
uA = function(t,alpha,beta,u0){
  return(alpha / beta * (1 - exp(- beta * t)) + u0 * exp(- beta * t))
}
# analytical solution for spliced mRNA
sA = function(t,alpha,beta,gamma,u0,s0){
  egt = exp(-gamma*t)
  if(beta-gamma==0){
    return(s0 * egt + alpha/gamma * (1-egt) + (beta * u0-alpha)*t*egt)
  }
  return(s0 * egt + alpha/gamma * (1-egt) + (beta * u0-alpha)/(beta - gamma)*(egt-exp(-beta * t)))
}
# calculate correlation between transcription activity and mRNA abundance
trans_act = function(t,t1,t2,t3){
  if(t<t1|t>t2){
    return(0)
  }else{
    return(1)
  }
}

# Type 1: Up
# t: 0-t1, alpha = alpha1; t: t1-t2, alpha = alpha2
Type_1_fun = function(t1,t2,t3,deltat = 0.01,alpha1,alpha2,beta,gamma,u0= alpha1/beta,s0= alpha1/gamma){
  # part 1 t:0-t1, alpha = alpha1
  t.part1 = seq(0,t1-deltat,deltat)
  u.part1 = uA(t.part1,alpha = alpha1,beta,u0)
  s.part1 = sA(t.part1,alpha = alpha1,beta,gamma,u0,s0)
  
  # part 2 t:t1-t2, alpha = alpha2
  u0.part2 = u.part1[length(u.part1)]
  s0.part2 = s.part1[length(s.part1)]
  t.part2 = seq(t1,t2-deltat,deltat)
  u.part2 = uA(t.part2-t1,alpha = alpha2,beta,u0)
  s.part2 = sA(t.part2-t1,alpha = alpha2,beta,gamma,u0,s0)
  
  # combine 2 parts
  t = c(t.part1,t.part2)
  ut = c(u.part1,u.part2)
  st = c(s.part1,s.part2)
  
  # return dataframe
  df = data.frame(t = t, ut = ut, st = st)
  return(df)
} 
# Type 2: Down
# t: 0-t1, alpha = alpha2; t: t1-t2, alpha = alpha1
Type_2_fun = function(t1,t2,t3,deltat = 0.01,alpha1,alpha2,beta,gamma,u0= alpha2/beta,s0= alpha2/gamma){
  # part 1 t:0-t1, alpha = alpha1
  t.part1 = seq(0,t1-deltat,deltat)
  u.part1 = uA(t.part1,alpha = alpha2,beta,u0)
  s.part1 = sA(t.part1,alpha = alpha2,beta,gamma,u0,s0)
  
  # part 2 t:t1-t2, alpha = alpha2
  u0.part2 = u.part1[length(u.part1)]
  s0.part2 = s.part1[length(s.part1)]
  t.part2 = seq(t1,t2-deltat,deltat)
  u.part2 = uA(t.part2-t1,alpha = alpha1,beta,u0)
  s.part2 = sA(t.part2-t1,alpha = alpha1,beta,gamma,u0,s0)
  
  # combine 2 parts
  t = c(t.part1,t.part2)
  ut = c(u.part1,u.part2)
  st = c(s.part1,s.part2)
  
  # return dataframe
  df = data.frame(t = t, ut = ut, st = st)
  return(df)
} 
# Type 3: Up and Down
# t: 0-t1, alpha = alpha1; t: t1-t2, alpha = alpha2; t: t2-t3, alpha = alpha1
Type_3_fun = function(t1,t2,t3,deltat = 0.01,alpha1,alpha2,beta,gamma,u0= alpha1/beta,s0= alpha1/gamma){
  # part 1 t:0-t1, alpha = alpha1
  t.part1 = seq(0,t1-deltat,deltat)
  u.part1 = uA(t.part1,alpha = alpha1,beta,u0)
  s.part1 = sA(t.part1,alpha = alpha1,beta,gamma,u0,s0)
  
  # part 2 t:t1-t2, alpha = alpha2
  u0.part2 = u.part1[length(u.part1)]
  s0.part2 = s.part1[length(s.part1)]
  t.part2 = seq(t1,t2-deltat,deltat)
  u.part2 = uA(t.part2-t1,alpha = alpha2,beta,u0)
  s.part2 = sA(t.part2-t1,alpha = alpha2,beta,gamma,u0,s0)
  
  # part 3 t: t2-t3, alpha = alpha1
  u0.part3 = u.part2[length(u.part2)]
  s0.part3 = s.part2[length(s.part2)]
  t.part3 = seq(t2,t3-deltat,deltat)
  u.part3 = uA(t.part3-t2,alpha = alpha1,beta,u0.part3)
  s.part3 = sA(t.part3-t2,alpha = alpha1,beta,gamma,u0.part3,s0.part3)
  
  # combine 3 parts
  t = c(t.part1,t.part2,t.part3)
  ut = c(u.part1,u.part2,u.part3)
  st = c(s.part1,s.part2,s.part3)
  
  # return dataframe
  df = data.frame(t = t, ut = ut, st = st)
  return(df)
}
# Type 4: Down and Up
# t: 0-t1, alpha = alpha2; t: t1-t2, alpha = alpha1; t: t2-t3, alpha = alpha2
Type_4_fun = function(t1,t2,t3,deltat = 0.01,alpha1,alpha2,beta,gamma,u0= alpha2/beta,s0= alpha2/gamma){
  # part 1 t:0-t1, alpha = alpha1
  t.part1 = seq(0,t1-deltat,deltat)
  u.part1 = uA(t.part1,alpha = alpha2,beta,u0)
  s.part1 = sA(t.part1,alpha = alpha2,beta,gamma,u0,s0)
  
  # part 2 t:t1-t2, alpha = alpha2
  u0.part2 = u.part1[length(u.part1)]
  s0.part2 = s.part1[length(s.part1)]
  t.part2 = seq(t1,t2-deltat,deltat)
  u.part2 = uA(t.part2-t1,alpha = alpha1,beta,u0)
  s.part2 = sA(t.part2-t1,alpha = alpha1,beta,gamma,u0,s0)
  
  # part 3 t: t2-t3, alpha = alpha1
  u0.part3 = u.part2[length(u.part2)]
  s0.part3 = s.part2[length(s.part2)]
  t.part3 = seq(t2,t3-deltat,deltat)
  u.part3 = uA(t.part3-t2,alpha = alpha2,beta,u0.part3)
  s.part3 = sA(t.part3-t2,alpha = alpha2,beta,gamma,u0.part3,s0.part3)
  
  # combine 3 parts
  t = c(t.part1,t.part2,t.part3)
  ut = c(u.part1,u.part2,u.part3)
  st = c(s.part1,s.part2,s.part3)
  
  # return dataframe
  df = data.frame(t = t, ut = ut, st = st)
  return(df)
}

# calculate accuracy
cal_accuracy = function(a,b){
  df = cbind(a,b)
  TP = length(which(rowSums(df)==2))
  TN = length(which(rowSums(df)==0))
  accuracy = (TP+TN)/nrow(df)
  return(accuracy)
}

# parameters
if(T){
  regtime = 1000
  t1 = regtime
  deltat2 = regtime
  deltat3 = regtime
  t2 = t1 + deltat2
  t3 = t2 + deltat3
  deltat = 0.01
  alpha1 = 1
  alpha2 = 10
  beta = 0.1
  gamma = 0.01
}
par(mfrow = c(2,2))
type_id = 3
for(type_id in 1:4){
  # calculate
  if(T){
    type_funs = list(Type_1_fun,Type_2_fun,Type_3_fun,Type_4_fun)
    result = type_funs[[type_id]](t1,t2,t3,deltat,alpha1,alpha2,beta,gamma)
    t = result$t
    ut = result$ut
    st = result$st
    trans_a = sapply(t,trans_act,t1,t2,t3)
    if(type_id ==2|type_id==4){trans_a = 1-trans_a}  # change active or inactive 
    
    # calculate accuracy
    thr = (alpha1+alpha2)/2/beta
    ut.bina = as.numeric(ut>=thr)
    thr = (alpha1+alpha2)/2/gamma
    st.bina = as.numeric(st>=thr)
    
    u.acc = cal_accuracy(ut.bina,trans_a)
    s.acc = cal_accuracy(st.bina,trans_a)
    # theoretical accuracy
    if(type_id ==1|type_id==2){
      u.acc.thr = 1- log(2)/beta/t2
      s.acc.thr = 1- log(2)/gamma/t2
    }else{
      u.acc.thr = 1- 2*log(2)/beta/t3
      s.acc.thr = 1- 2*log(2)/gamma/t3
    }
      
    
  }
  # plot
  pdf("mRNA dynamics ODE.pdf",width = 7, height = 7)
  if(T){
    # mRNA dynamics
    plot(t,ut,ylim = c(0,max(st,ut)),type = "l",col = "red",ylab = "mRNA",
         main = paste("mRNA dynamics: regulation time =",regtime,sep = " "))
    points(t,st,col = "blue",type = "l")
    points(t,max(st,ut)*trans_a,col = "purple",type = "l",lty = 2)
    cols = c("red","blue")
    
    # legend("topleft",c("u ( t )","s ( t )"),text.col = cols,col = cols,lty = 1)
    # abline(h = (alpha1+alpha2)/2/beta,col = "red",lty = 2)
    # abline(h = (alpha1+alpha2)/2/gamma,col = "blue",lty = 2)
    
  }
  dev.off()
  
  # plot accuracy
  if(F){
    barplot(height = c(u.acc,u.acc.thr,s.acc,s.acc.thr),names.arg = c("pre-mRNA","theory","mRNA","theory"),
            col = c("red","red","blue","blue"),ylim = c(0,1),main = "accuracy")
    
    abline(h = 1,col = "purple",lty = 2)
  }
}

# Use Type 3 function: Up and Down
# Investigate the effect of regulation time (regtime)
# consider varied timescale of transcription regulation
type_id = 3
par(mfrow = c(2,2))
type_funs = list(Type_1_fun,Type_2_fun,Type_3_fun,Type_4_fun)
for(type_id in 1:4){
  df = NULL
  for(regtime in 10**seq(2,4,0.25)){
    t1 = regtime
    deltat2 = regtime
    deltat3 = regtime
    t2 = t1 + deltat2
    t3 = t2 + deltat3
    deltat = 0.1
    result = type_funs[[type_id]](t1,t2,t3,deltat,alpha1,alpha2,beta,gamma)
    t = result$t
    ut = result$ut
    st = result$st
    trans_a = sapply(t,trans_act,t1,t2,t3)
    if(type_id ==2|type_id==4){trans_a = 1-trans_a}  # change active or inactive 
    
    # calculate accuracy
    thr = (alpha1+alpha2)/2/beta
    ut.bina = as.numeric(ut>=thr)
    thr = (alpha1+alpha2)/2/gamma
    st.bina = as.numeric(st>=thr)
    
    u.acc = cal_accuracy(ut.bina,trans_a)
    s.acc = cal_accuracy(st.bina,trans_a)
    # theoretical accuracy
    if(type_id ==1|type_id==2){
      u.acc.thr = 1- log(2)/beta/t2
      s.acc.thr = 1- log(2)/gamma/t2
    }else{
      u.acc.thr = 1- 2*log(2)/beta/t3
      s.acc.thr = 1- 2*log(2)/gamma/t3
    }
    
    df = rbind(df,c(regtime,u.acc,s.acc,u.acc.thr,s.acc.thr))
    print(paste(regtime,sep = " "))
    
    # plot mRNA dynamics
    if(T){
      # mRNA dynamics
      plot(t,ut,ylim = c(0,max(st,ut)),type = "l",col = "red",ylab = "mRNA",
           main = paste("mRNA dynamics: regulation time =",regtime,sep = " "))
      points(t,st,col = "blue",type = "l")
      points(t,max(st,ut)*trans_a,col = "purple",type = "l",lty = 2)
      cols = c("red","blue")
      
      # legend("topleft",c("u ( t )","s ( t )"),text.col = cols,col = cols,lty = 1)
      abline(h = (alpha1+alpha2)/2/beta,col = "red",lty = 2)
      abline(h = (alpha1+alpha2)/2/gamma,col = "blue",lty = 2)
      
    }
  }
  df = as.data.frame(df)
  colnames(df) = c("regtime","u.acc","s.acc","u.acc.thr","s.acc.thr")
  
  
  pdf("accuracy_regtime.pdf",width = 7, height = 7)
  
  # plot for df accuracy
  if(T){
    par(mfrow = c(1,1))
    t_alpha = df[,1]
    plot(t_alpha,df$u.acc,col = "red", type = "b",main = "effect of regulation time", log = "x",ylim = c(0.6,1),
         xlab = "regulation time", ylab = "accuracy")
    points(t_alpha,df$s.acc,col = "blue",type = "b")
    abline(h = 1,col = "purple",lty = 2)
    
    # points(t_alpha,df$u.acc.thr,col = "red", type = "b",ylim = c(0,1),pch = 2)
    # points(t_alpha,df$s.acc.thr,col = "blue",type = "b",pch = 2)
    
    # cols = c("red","red","blue","blue")
    # pchs = c(1,2,1,2)
    # legend("bottomright",c("u","u.thr","s","s.thr"),text.col = cols,col = cols,lty = 1,pch = pchs)
    cols = c("red","blue")
    legend("bottomright",c("u","s"),text.col = cols,col = cols,lty = c(1,1),inset = 0.05,cex = 1.3)
  }
  dev.off()
}

# Investigate the effect of parameters
# parameters: beta, gamma: 2**0*0.01,...,2**10*0.01
df = NULL
for(beta in 2**seq(-10,1,1)){
  for(gamma in 2**seq(-10,1,1)){
    type_id = 3
    result = type_funs[[type_id]](t1,t2,t3,deltat,alpha1,alpha2,beta,gamma)
    t = result$t
    ut = result$ut
    st = result$st
    trans_a = sapply(t,trans_act,t1,t2,t3)
    if(type_id ==2|type_id==4){trans_a = 1-trans_a}  # change active or inactive 
    
    # calculate accuracy
    thr = (alpha1+alpha2)/2/beta
    ut.bina = as.numeric(ut>=thr)
    thr = (alpha1+alpha2)/2/gamma
    st.bina = as.numeric(st>=thr)
    
    u.acc = cal_accuracy(ut.bina,trans_a)
    s.acc = cal_accuracy(st.bina,trans_a)
    # theoretical accuracy
    if(type_id ==1|type_id==2){
      u.acc.thr = 1- log(2)/beta/t2
      s.acc.thr = 1- log(2)/gamma/t2
    }else{
      u.acc.thr = 1- 2*log(2)/beta/t3
      s.acc.thr = 1- 2*log(2)/gamma/t3
    }
    
    df = rbind(df,c(beta,gamma,u.acc,s.acc,u.acc.thr,s.acc.thr))
    print(paste(beta,gamma,sep = " "))
  }
}
df = as.data.frame(df)
colnames(df) = c("beta","gamma","u.acc","s.acc","u.acc.thr","s.acc.thr")
# ggplot for accuracy
if(T){
  library(ggplot2)
  u.acc = df$u.acc
  s.acc = df$s.acc
  
  # plot for u.cor and s.cor
  project = log10(df[,1:2])
  colnames(project) = c("log10_beta","log10_gamma")
  project.new <- cbind(as.data.frame(project), u.acc=u.acc,s.acc = s.acc)
  head(project.new)
  # plot for u.acc
  p = ggplot(project.new, aes(x=log10_gamma, y=log10_beta,color = u.acc)) + geom_point(size = 10)
  p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
  p
  # plot for s.acc
  p = ggplot(project.new, aes(x=log10_gamma, y=log10_beta,color = s.acc)) + geom_point(size = 10)
  p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
  p
  # plot for u.acc - s.acc
  dif = u.acc - s.acc
  p = ggplot(project.new, aes(x=log10_gamma, y=log10_beta,color = dif)) + geom_point(size = 10)
  p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
  p
  
  ###  compare simulation results with theory
  # Effect of beta
  plot(df$beta,df$s.acc,col = "blue",log = "x",ylim = c(0,1),xlab = "beta",ylab = "accuracy",main ="Effect of beta")
  # lines(df$beta,1 - 2*log(2)/df$beta/t3,col = "red",lty = 2,pch = 2,type = "b")
  points(df$beta,df$u.acc,col = "red",type = "b")
  abline(h = 1,col = "purple",lty = 2)
  
  # Effect of gamma
  plot(df$gamma,df$s.acc,col = "blue",log = "x",ylim = c(0,1),xlab = "gamma",ylab = "accuracy",main ="Effect of gamma")
  lines(unique(df$gamma),1 - 2*log(2)/unique(df$gamma)/t3,col = "blue",lty = 2,pch = 2,type = "b")
  points(df$gamma,df$u.acc,col = "red",type = "p")
  abline(h = 1,col = "purple",lty = 2)
  
  # Select subset
  id = which(df$beta>2*df$gamma)
  plot(df$beta[id],df$s.acc[id],col = "blue",log = "x",ylim = c(0,1),xlab = "beta",ylab = "accuracy",main ="Effect of beta")
  points(df$beta[id],df$s.acc.thr[id],col = "blue",pch = 2)
  
}


