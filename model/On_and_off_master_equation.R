# On and Off model -- master equation simulaiton   # 2021/5/1
library(GillespieSSA)
library(infotheo)
# calculate correlation between transcription activity and mRNA abundance
trans_act = function(t,t1,t2,t3){
  if(t<t1|t>t2){
    return(0)
  }else{
    return(1)
  }
}
# SSA: Type 3: Up and Down
# t: 0-t1, alpha = alpha1; t: t1-t2, alpha = alpha2; t: t2-t3, alpha = alpha1
On_and_off_SSA = function(t1,t2,t3,deltat = 1,alpha1,alpha2,beta,gamma,u0= alpha1/beta,s0= alpha1/gamma,a,nu,seed_id){
  set.seed(seed_id)
  # part 1 t:0-t1, alpha = alpha1
  parms <- c(alpha = alpha1, beta = beta, gamma = gamma)
  x0 <- c(X1 = u0, X2 = s0,X3 = 1) 
  out <- ssa(x0 = x0,a = a,nu = nu,parms = parms,tf = t1,method = ssa.etl(tau = deltat),verbose = F,consoleInterval = 1) 
  t.part1 = seq(0,t1-deltat,deltat)
  u.part1 = out$data[,2][1:length(t.part1)]
  s.part1 = out$data[,3][1:length(t.part1)]
  u0 = out$data[,2][length(t.part1)+1]
  s0 = out$data[,3][length(t.part1)+1]
  
  # part 2 t:t1-t2, alpha = alpha2
  parms <- c(alpha = alpha2, beta = beta, gamma = gamma)
  x0 <- c(X1 = u0, X2 = s0,X3 = 1) 
  out <- ssa(x0 = x0,a = a,nu = nu,parms = parms,tf = t2-t1,method = ssa.etl(tau = deltat),verbose = F,consoleInterval = 1) 
  t.part2 = seq(t1,t2-deltat,deltat)
  u.part2 = out$data[,2][1:length(t.part2)]
  s.part2 = out$data[,3][1:length(t.part2)]
  u0 = out$data[,2][length(t.part2)+1]
  s0 = out$data[,3][length(t.part2)+1]
  
  # part 3 t: t2-t3, alpha = alpha1
  parms <- c(alpha = alpha1, beta = beta, gamma = gamma)
  x0 <- c(X1 = u0, X2 = s0,X3 = 1) 
  out <- ssa(x0 = x0,a = a,nu = nu,parms = parms,tf = t3-t2,method = ssa.etl(tau = deltat),verbose = F,consoleInterval = 1) 
  t.part3 = seq(t2,t3-deltat,deltat)
  u.part3 = out$data[,2][1:length(t.part3)]
  s.part3 = out$data[,3][1:length(t.part3)]
  
  # combine 3 parts
  t = c(t.part1,t.part2,t.part3)
  ut = c(u.part1,u.part2,u.part3)
  st = c(s.part1,s.part2,s.part3)
  
  # return dataframe
  df = data.frame(t = t, ut = ut, st = st)
  return(df)
}
# SSA ordinary
On_and_off_SSA_d = function(t1,t2,t3,deltat = 1,alpha1,alpha2,beta,gamma,u0= alpha1/beta,s0= alpha1/gamma,a,nu,seed_id){
  set.seed(seed_id)
  # part 1 t:0-t1, alpha = alpha1
  parms <- c(alpha = alpha1, beta = beta, gamma = gamma)
  x0 <- c(X1 = u0, X2 = s0,X3 = 1) 
  out <- ssa(x0 = x0,a = a,nu = nu,parms = parms,tf = t1,method = ssa.d(),verbose = F,consoleInterval = 1) 
  # change to discrete time
  result = change_time(out$data,deltat,t1)
  t.part1 = seq(0,t1-deltat,deltat)
  u.part1 = result[,2][1:length(t.part1)]
  s.part1 = result[,3][1:length(t.part1)]
  u0 = result[,2][length(t.part1)]
  s0 = result[,3][length(t.part1)]
  
  # part 2 t:t1-t2, alpha = alpha2
  parms <- c(alpha = alpha2, beta = beta, gamma = gamma)
  x0 <- c(X1 = u0, X2 = s0,X3 = 1) 
  out <- ssa(x0 = x0,a = a,nu = nu,parms = parms,tf = t2-t1,method = ssa.d(),verbose = F,consoleInterval = 1)  
  result = change_time(out$data,deltat,t2-t1)
  t.part2 = seq(t1,t2-deltat,deltat)
  u.part2 = result[,2][1:length(t.part2)]
  s.part2 = result[,3][1:length(t.part2)]
  u0 = result[,2][length(t.part2)]
  s0 = result[,3][length(t.part2)]
  
  # part 3 t: t2-t3, alpha = alpha1
  parms <- c(alpha = alpha1, beta = beta, gamma = gamma)
  x0 <- c(X1 = u0, X2 = s0,X3 = 1) 
  out <- ssa(x0 = x0,a = a,nu = nu,parms = parms,tf = t3-t2,method = ssa.d(),verbose = F,consoleInterval = 1) 
  result = change_time(out$data,deltat,t3-t2)
  t.part3 = seq(t2,t3-deltat,deltat)
  u.part3 = result[,2][1:length(t.part3)]
  s.part3 = result[,3][1:length(t.part3)]
  
  # combine 3 parts
  t = c(t.part1,t.part2,t.part3)
  ut = c(u.part1,u.part2,u.part3)
  st = c(s.part1,s.part2,s.part3)
  
  # return dataframe
  df = data.frame(t = t, ut = ut, st = st)
  return(df)
}
change_time = function(df, deltat,tmax){
  result = matrix(0,nrow = tmax/deltat,ncol = ncol(df))
  tpoints = df[,1]
  for(i in 1:length(tpoints)){
    tstart = floor(tpoints[i]/deltat)+1
    if(i!=length(tpoints)){
      tend = floor(tpoints[i+1]/deltat)+1
    }else{
      tend = tmax/deltat
    }
    if(tend>tmax/deltat){ 
      tend = tmax/deltat
      result[tstart:tend,] = matrix(rep(df[i,],length(tstart:tend)),nrow= length(tstart:tend),byrow = T)
      return(result)
    }
    result[tstart:tend,] = matrix(rep(df[i,],length(tstart:tend)),nrow= length(tstart:tend),byrow = T)
  }
  return(result)
}
# ODE: Type 3: Up and Down
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
  deltat = 0.1
  
  # Define propensity functions
  a  <- c("alpha", "beta*X1","gamma*X2")
  # Define state-change matrix
  nu <- matrix(c(+1, -1, 0, 0, 1, -1,0,0,0), nrow = 3, byrow = TRUE)
  alpha2 = 0.1
  alpha1 = alpha2/10
  beta = 0.1
  gamma = 0.01
  u0= ceiling(alpha1/beta)
  s0= ceiling(alpha1/gamma)
}
alphas = 10**seq(-1,1,0.25)
alphas = 10**c(seq(-1,-0.1,0.1),seq(0,1,0.25))

df = NULL
for(alpha_id in 1:length(alphas)){
  alpha = alphas[alpha_id]
  alpha2 = alpha
  alpha1 = alpha2/10
  u0= ceiling(alpha1/beta)
  s0= ceiling(alpha1/gamma)
  # calculate
  # seed_id = 1
  type_id = 3
  df.tmp = NULL
  for(seed_id in 1:10){
    set.seed(seed_id)
    if(T){
      result = On_and_off_SSA_d(t1,t2,t3,deltat,alpha1,alpha2,beta,gamma,u0,s0,a,nu,seed_id)
      t = result$t
      ut = result$ut
      st = result$st
      trans_a = sapply(t,trans_act,t1,t2,t3)
      
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
    # calculate mutual information
    if(T){
      X = data.frame(x = trans_a, y = ut, z = st)
      mi.all = mutinformation(X)
      mi.gg = mi.all[1,1]
      mi.gu = mi.all[1,2]
      mi.gs = mi.all[1,3]
    }
    df.tmp = rbind(df.tmp,c(u.acc,s.acc,mi.gg,mi.gu,mi.gs))
    
  }
  df = rbind(df,apply(df.tmp,2,mean))
  print(alpha_id)
}
colnames(df) = c("u.acc","s.acc","mi.gg","mi.gu","mi.gs")
df = data.frame(df)

# plot
# pdf("mRNA dynamics stochastic.pdf",width = 7, height = 7)
if(T){
  # mRNA dynamics
  plot(t,ut,ylim = c(0,max(st,ut)),type = "l",col = "red",ylab = "mRNA",
       main = paste("mRNA dynamics: alpha =",alpha2,sep = " "))
  points(t,st,col = "blue",type = "l")
  points(t,max(st,ut)*trans_a,col = "purple",type = "l",lty = 2)
  cols = c("red","blue")
  
  # legend("topleft",c("u ( t )","s ( t )"),text.col = cols,col = cols,lty = 1)
  # abline(h = (alpha1+alpha2)/2/beta,col = "red",lty = 2)
  # abline(h = (alpha1+alpha2)/2/gamma,col = "blue",lty = 2)
  
}
# dev.off()
# plot MI(t): t1<t<t2
if(T){
  ut.all = NULL
  st.all = NULL
  for(seed_id in 1:10){
    set.seed(seed_id)
    if(T){
      result = On_and_off_SSA(t1,t2,t3,deltat,alpha1,alpha2,beta,gamma,u0,s0,a,nu,seed_id)
      t = result$t
      ut = result$ut
      st = result$st
      trans_a = sapply(t,trans_act,t1,t2,t3)
    }
    ut.all = rbind(ut.all,ut)
    st.all = rbind(st.all,st)
  }
  timepoint = 3
  X = data.frame(x = rep(trans_a[timepoint],nrow(ut.all)), y = ut.all[,timepoint], z = st.all[,timepoint])
  mutinformation(X)
  mutinformation(X)[1,2]
  mutinformation(X)[1,3]
  samplepoint = 2
  X = data.frame(x = trans_a, y = ut.all[samplepoint,], z = st.all[samplepoint,])
  mutinformation(X)
  mi.gg = mutinformation(X)[1,1]
  mi.gu = mutinformation(X)[1,2]
  mi.gs = mutinformation(X)[1,3]
  
}

# plot for determined model
if(T){
  result = Type_3_fun(t1,t2,t3,deltat,alpha1,alpha2,beta,gamma)
  t = result$t
  ut = result$ut
  st = result$st
  points(t,ut,col = "red",type = "l")
  points(t,st,col = "blue",type = "l")
  
  # calculate accuracy
  thr = (alpha1+alpha2)/2/beta
  ut.bina = as.numeric(ut>=thr)
  thr = (alpha1+alpha2)/2/gamma
  st.bina = as.numeric(st>=thr)
  u.acc.ode = cal_accuracy(ut.bina,trans_a)
  s.acc.ode = cal_accuracy(st.bina,trans_a)
  
}
# plot accuracy
if(T){
  barplot(height = c(u.acc,u.acc.ode,u.acc.thr,s.acc,s.acc.ode,s.acc.thr),
          names.arg = c("pre-mRNA","ode","theory","mRNA","ode","theory"),
          col = c("red","red","red","blue","blue","blue"),ylim = c(0,1),main = "accuracy")
  
  abline(h = 1,col = "purple",lty = 2)
}

pdf("accuracy_alpha.pdf",width = 7, height = 7)
# plot accuracty for different alphas
if(T){
  df.tmp = cbind(df,alphas)
  plot(alphas,df$u.acc, col = "red",ylim = c(0.8,1),log = 'x',type= "b",xlab = "alpha",ylab = "accuracy",
       main = "effect of transcription rate")
  points(alphas,df$s.acc, col = "blue",ylim = c(0,1),type= "b")
  # lines(alphas,predict(loess(df$u.acc~alphas)),col = "red")
  # abline(lm(s.acc~alphas,data = df.tmp),col= "blue")
  abline(h = 1,col = "purple",lty = 2)
  cols = c("red","blue")
  legend("bottomright",c("u","s"),text.col = cols,col = cols,lty = c(1,1),inset = 0.05,cex = 1.3)
  
}
dev.off()
# plot MI for different alphas
if(T){
  plot(alphas,df$mi.gu, col = "red",ylim = c(0,1),log = 'x',type= "b",xlab = "alpha",ylab = "MI")
  points(alphas,df$mi.gs, col = "blue",ylim = c(0,1),type= "b")
  abline(h = df$mi.gg,col = "purple")
}
# Calculate MI for different sample number
df = NULL
if(T){
  deltat = 1
  alpha = 1
  alpha2 = alpha
  alpha1 = alpha2/10
  u0= ceiling(alpha1/beta)
  s0= ceiling(alpha1/gamma)
  # calculate
  # seed_id = 1
  type_id = 3
  u.all = NULL
  s.all = NULL
  df.tmp = NULL
  for(seed_id in 1:10){
    set.seed(seed_id)
    if(T){
      result = On_and_off_SSA_d(t1,t2,t3,deltat,alpha1,alpha2,beta,gamma,u0,s0,a,nu,seed_id)
      t = result$t
      ut = result$ut
      st = result$st
      trans_a = sapply(t,trans_act,t1,t2,t3)
      
      u.all = rbind(u.all,ut)
      s.all = rbind(s.all,st)
      
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
    # calculate mutual information
    if(T){
      X = data.frame(x = rep(trans_a,seed_id), y = as.vector(t(u.all)), z = as.vector(t(st)))
      mi.all = mutinformation(X)
      mi.gg = mi.all[1,1]
      mi.gu = mi.all[1,2]
      mi.gs = mi.all[1,3]
    }
    df.tmp = rbind(df.tmp,c(u.acc,s.acc,mi.gg,mi.gu,mi.gs))
    print(seed_id)
  }
  colnames(df.tmp) = c("u.acc","s.acc","mi.gg","mi.gu","mi.gs")
  df.tmp = data.frame(df.tmp)
  
  plot(1:seed_id,df.tmp$mi.gu,col = "red")
  
}
# plot dynamics
if(T){
  par(mfrow =c(1,2))
  plot(t,u.all[1,],type = "l",col= "red")
  lines(t,u.all[2,],type = "l",col= "red")
  lines(t,u.all[3,],type = "l",col= "red")
  lines(t,u.all[4,],type = "l",col= "red")
  lines(t,u.all[5,],type = "l",col= "red")
  
  plot(t,s.all[1,],type = "l",col= "blue")
  lines(t,s.all[2,],type = "l",col= "blue")
  lines(t,s.all[3,],type = "l",col= "blue")
  lines(t,s.all[4,],type = "l",col= "blue")
  lines(t,s.all[5,],type = "l",col= "blue")
  
  
  hist(as.vector(u.all))
  hist(as.vector(s.all))
  
}


# effect of both alpha and regulation time
df = NULL
alphas = 10**seq(-1,0,0.1)
regtimes = 10**seq(2,3,0.1)
for(alpha in alphas){
  alpha2 = alpha
  alpha1 = alpha2/10
  u0= ceiling(alpha1/beta)
  s0= ceiling(alpha1/gamma)
  # calculate
  # seed_id = 1
  type_id = 3
  for(regtime in regtimes){
    t1 = regtime
    deltat2 = regtime
    deltat3 = regtime
    t2 = t1 + deltat2
    t3 = t2 + deltat3
    deltat = 0.1
    
    type_id = 3
    df.tmp = NULL
    for(seed_id in 1:5){
      set.seed(seed_id+10)
      if(T){
        result = On_and_off_SSA_d(t1,t2,t3,deltat,alpha1,alpha2,beta,gamma,u0,s0,a,nu,seed_id)
        t = result$t
        ut = result$ut
        st = result$st
        trans_a = sapply(t,trans_act,t1,t2,t3)
        
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
      # calculate mutual information
      if(T){
        X = data.frame(x = trans_a, y = ut, z = st)
        mi.all = mutinformation(X)
        mi.gg = mi.all[1,1]
        mi.gu = mi.all[1,2]
        mi.gs = mi.all[1,3]
      }
      df.tmp = rbind(df.tmp,c(alpha, regtime,u.acc,s.acc,mi.gg,mi.gu,mi.gs))
    }
    
    df = rbind(df,apply(df.tmp,2,mean))
    print(paste(alpha, regtime))
  }
}
df = as.data.frame(df)
colnames(df) = c("alpha", "regtime", "u.acc","s.acc","mi.gg","mi.gu","mi.gs")

# ggplot for accuracy
if(T){
  library(ggplot2)
  u.acc = df$u.acc
  s.acc = df$s.acc
  
  u.acc = df$mi.gu
  s.acc = df$mi.gs
  
  # plot for u.cor and s.cor
  project = log10(df[,1:2])
  colnames(project) = c("log10_alpha","log10_regtime")
  project.new <- cbind(as.data.frame(project), u.acc=u.acc,s.acc = s.acc)
  head(project.new)
  # plot for u.acc
  p = ggplot(project.new, aes(x=log10_alpha, y=log10_regtime,color = u.acc)) + geom_point(size = 10)
  p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
  p
  # plot for s.acc
  p = ggplot(project.new, aes(x=log10_alpha, y=log10_regtime,color = s.acc)) + geom_point(size = 10)
  p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
  p
  # plot for u.acc - s.acc
  dif = u.acc - s.acc
  p = ggplot(project.new, aes(x=log10_alpha, y=log10_regtime,color = dif)) + geom_point(size = 10)
  p = p + scale_color_gradient2(low = "blue", mid = "white", high = "red")
  p
  
  ###  compare simulation results with theory
  # Effect of alpha
  plot(df$alpha,df$s.acc,col = "blue",log = "x",ylim = c(0,1),xlab = "alpha",ylab = "accuracy",main ="Effect of alpha")
  points(df$alpha,df$u.acc,col = "red")
  abline(h = 1,col = "purple",lty = 2)
  
  # Effect of regulation time
  plot(df$regtime,df$s.acc,col = "blue",log = "x",ylim = c(0,1),xlab = "regtime",ylab = "accuracy",main ="Effect of regtime")
  points(df$regtime,df$u.acc,col = "red",type = "p")
  abline(h = 1,col = "purple",lty = 2)
  
  ###  compare simulation results with theory
  # Effect of alpha
  plot(df$alpha,df$mi.gs,col = "blue",log = "x",ylim = c(0,0.8),xlab = "alpha",ylab = "MI",main ="Effect of alpha")
  points(df$alpha,df$mi.gu,col = "red")
  abline(h = df$mi.gg,col = "purple",lty = 2)
  
  # Effect of regulation time
  plot(df$regtime,df$mi.gs,col = "blue",log = "x",ylim = c(0,0.8),xlab = "regtime",ylab = "MI",main ="Effect of regtime")
  points(df$regtime,df$mi.gu,col = "red",type = "p")
  abline(h = df$mi.gg,col = "purple",lty = 2)
  
  
}
