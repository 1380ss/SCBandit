#Goal: givne output, 
#load that function TS2b
#1. function calculating optimization of objective with differnt lambda
# look at your overleaf..
#2. calculate empirical performance (at each step, each 'c' and epsilon.)
library(progress)
library(utils)
#https://r-coder.com/progress-bar-r/ 

library(ggplot2)
library(patchwork)
source("Algorithm for two arm.R")
source("Algorithms for general use.R")


################ Algorithm for simulation and plot

#functions needed
#function calculating optimal solution
entropy_objective_optimal_2arm <- function(lambda,delta,x=NULL,min_sacle=0){
  #delta = mu1-mu2, x is proportion to arm 1
  #objective = delta x + lambda [xlogx + (1-x) log(1-x)]
  xx <- exp(delta/lambda)
  if(is.null(x)){
    x <- xx/(1+xx)
  }
  scale_factor <- (lambda*log(1+xx))
  if(is.na(min_sacle)){
    min_sacle <- delta*0.5+lambda*log(2)
  }
  
  loss <- (x*delta-lambda*( x*log(x)+(1-x)*log(1-x))-min_sacle)/(scale_factor-min_sacle)
  return(list(best_x=xx/(1+xx),loss=loss))
}
# entropy_objective_optimal_2arm(lambda=0.2,delta=0.2,x=0.9)
create_df_setting <- function(parameter_list,effect_list){
  algorithm_list <- names(parameter_list)
  df1 <- data.frame()
  for(i in algorithm_list){
    df1 <- rbind(df1,expand.grid(algorithm=as.character(i),parameter=parameter_list[[i]]))
  }
  df1$algorithm <- as.character(df1$algorithm )
  df2 <- data.frame(effect=effect_list)#at this point you only have effect size. 
  
  df_settings <- merge(df1,df2)
  df_settings$algo_name <- df_settings$algorithm #name with parameter, so can be used for plots
  df_settings$algo_name[!is.na(df_settings$parameter)] <- paste0(df_settings$algo_name[!is.na(df_settings$parameter)],' (', df_settings$parameter[!is.na(df_settings$parameter)],')')
  
  #process epsilon and c
  #can write a better work flow later
  df_settings$epsilon <- 0
  df_settings$c <- 0
  df_settings$epsilon[df_settings$algorithm=='UR'] <- 1
  df_settings$epsilon[df_settings$algorithm=='TT-TS'] <- df_settings$parameter[df_settings$algorithm=='TT-TS']
  
  df_settings$c[df_settings$algorithm=='UR'] <- 1
  df_settings$c[df_settings$algorithm=='TS PostDiff'] <- df_settings$parameter[df_settings$algorithm=='TS PostDiff']
  
  return(df_settings) 
}




############## main

#Setting <- 'pa=0.6,pb=0.4,n=785,c=0.1,eps=0'

effect_list <- c(0,0.05,0.1,0.15,0.2,0.3,0.4)
horizon_list <- c(100,200,500,785)
lambdas <- c(0.01,0.02,0.04,0.08,0.12,0.16,0.2)
#horizon_max <- max(horizon_list)

parameter_list <- list('TS'=c(NA),
                       'TS PostDiff'=c(0.05,0.075,0.1,0.15,0.2),
                       'TT-TS'=c(0.125,0.325),
                       'UR'=c(NA))
parameter_list <- list('TS'=c(),
                       'TS PostDiff'=c(0.05,0.075,0.1,0.15,0.2),
                       'TT-TS'=c(),
                       'UR'=c())

df_settings <- create_df_setting(parameter_list,effect_list)


B <- 10000
df_res <-df_settings
df_x <- array(dim=c(nrow(df_settings),length(horizon_list),B))
for(i in 1:nrow(df_settings)){
  s <- df_settings[i,]
  para <- list(pa=0.5+s$effect/2,pb=0.5-s$effect/2,n=max(horizon_list),eps=s$epsilon,c=s$c)
  res <- sim2(para,B=B,lambdas=lambdas)
  df_x[i,,] <- t(res$prop_x1[,horizon_list])
  #df_x[i,,] <- res$objective_scores[horizon_list,]
}


cp <- function(lambda,time_i=4,y='score'){
  df <- df_settings
  df$score <- 0
  if(y=='score'){
    for(i in 1:nrow(df_settings)){
      df$score[i] <- mean(entropy_objective_optimal_2arm (lambda,delta=df_settings[i,]$effect,x=df_x[i,time_i,])$loss)
      #can also plot x v.s. best x
    }
    summarize(group_by(df,algo_name),score_m=mean(score))
    ggplot(data=df,aes(x=effect,y=score,group=parameter,color=algo_name))+
      geom_line()+
      ylim(min=0.6,max=1)+
      ggtitle(paste0('lambda=',lambda,'; sample size =',horizon_list[time_i]))
  } else {
    for(i in 1:nrow(df_settings)){
      df$score[i] <- mean(x=df_x[i,time_i,])
      #can also plot x v.s. best x
    }
    summarize(group_by(df,algo_name),score_m=mean(score))
    ggplot(data=df,aes(x=effect,y=score,group=parameter,color=algo_name))+
      geom_line()+
      ylim(min=0.45,max=1)+
      ylab('allocation to arm1')+
      ggtitle(paste0('sample size =',horizon_list[time_i]))
    
  }
  
  
}
ccp(0.05)
lambdas <- c(0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2)
df <- df_settings

df$prop_x <- apply(df_x[,4,],1,mean)

df[df$parameter==0.15,][c(1,4,7),]
df$score <- 0
time_i <- 4
df0 <- summarize(group_by(df,algo_name),x=mean(score))[,1]
for (j in lambdas){
  for(i in 1:nrow(df_settings)){
    df$score[i] <- mean(entropy_objective_optimal_2arm (lambda=j,delta=df_settings[i,]$effect,x=df_x[i,time_i,])$loss)
    #can also plot x v.s. best x
  }
  df0 <- cbind(df0,summarize(group_by(df,algo_name),x=mean(score))[,2])
}
colnames(df0) <- c('algo',lambdas)
df0 <- rbind(df0,apply(df0,2,max))
write.csv(df0,file='TTTS.csv')
y <- c(0.025,0.05,0.075,0.075,0.1,0.1,0.15,0.15,0.2,0.2,0.2,0.25)
plot(lambdas,y)
abline(lm(y~lambdas))
summary(lm(y~lambdas))

y_TTTS <- c(0.025,0.025,0.05,0.1,0.1,0.2,0.3,0.3,0.3,0.5,0.6,0.7)
plot(lambdas,y_TTTS)
abline(lm(y_TTTS~lambdas))
summary(lm(y_TTTS~lambdas))
ccp <- function(lambda,y='score',time_list=c(2,4,7,9)){
  for (i in 1:4){
    assign(paste0('p',i),cp(lambda,time_i=time_list[i],y=y))
  }
  combined <- p1 + p2+p3+p4 & theme(legend.position = "bottom")
  combined + plot_layout(guides = "collect")
}

ccp(lambda=0.05)


# reminder for this function : df2table(df_res$algo_name,df_res$effect,df_x[,t,])

lambda <- 0.2
#x <- (1:100)/100
y <- exp(x/lambda)/(1+exp(x/lambda))



b <- 6

x <- (1:60)/100
mp <- function(b){
  cl <- c(0.075,0.1,0.15,0.3)
  colors <- c('green','blue','red','orange')
  c <- 0.025
  a <- (1+c*b-2*c)/(1-c)
  y <- pbeta(x,a,b)/2+0.5
  plot(x,y,main=paste0('b=',b))
  for(i in 1:length(cl)){
    c <- cl[i]
    
    a <- (1+c*b-2*c)/(1-c)
    y <- pbeta(x,a,b)/2+0.5
    lines(x,y,col=colors[i])
  }
  abline(v=0.3,h=0.7)
  abline(v=0.15)
}
mp(4)


## sim 2, fix lambda, differ effect size
plot_sim2 <- function(lambda,os_df,c,eps){
  lambda_index <- which(lambdas==lambda)
  df1 <- data.frame(score=c(t(os_df[,,lambda_index])),time_step=rep(c(1:para$n),length(ef)),effect=rep(as.character(ef),each=para$n))
  algo <- 'TS PostDiff with c = '
  algo_para <- c
  if(eps+c>=1){
    algo_para <- ''
    algo='UR'
  } else if(eps+c==0){
    algo_para <- ''
    algo='TS'
  } else if(eps>0){
    algo_para <- eps
    algo='TT-TS with epsilon = '
  }
  tit <- paste0(algo,algo_para,'; lambda=',lambda)
  ggplot(data=df1, aes(x=time_step, y=score, group=effect, color=effect)) +
    geom_line()+
    ggtitle(tit)
    #ylim(0.7,1)
}

#start
ef <- c(0,0.05,0.1,0.2,0.3)
lambdas <- c(0.05,0.1,0.2,0.3,0.5,1,2,3)
n <- 1000
c <- 0
eps <- 1
B <- 500
os_df <- array(dim=c(length(ef),n,length(lambdas)))
for (j in 1:length(ef)){
  para <- list(pa=0.5+ef[j]/2,pb=0.5-ef[j]/2,n=n,c=c,eps=eps)
  res <- sim2(para,B=B,lambdas=lambdas)
  os_df[j,,] <- res$objective_scores
}
p1 <- plot_sim2(lambda=0.05,os_df,c=0,eps=1)
p2 <- plot_sim2(lambda=0.1,os_df,c=0,eps=1)
p3 <- plot_sim2(lambda=0.2,os_df,c=0,eps=1)
p4 <- plot_sim2(lambda=0.3,os_df,c=0,eps=1)
p5 <- plot_sim2(lambda=0.5,os_df,c=0,eps=1)
p6 <- plot_sim2(lambda=1,os_df,c=0,eps=1)
p7 <- plot_sim2(lambda=2,os_df,c=0,eps=1)
p8 <- plot_sim2(lambda=3,os_df,c=0,eps=1)
combined <- p1 + p2+p3+p4+p5+p6+p7+p8 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
#lambdas <- c(0.05,0.1,0.2,0.3,0.5,1,2,3)
lambda <- 0.05
p1 <- plot_sim2(lambda,os_df_TS,c=0,eps=0)
p2 <- plot_sim2(lambda,os_df_PostDiff_01,c=0.1,eps=0)
p3 <- plot_sim2(lambda,os_df_TTTS_01,c=0,eps=0.1)

lambda <- 0.1
p4 <- plot_sim2(lambda,os_df_TS,c=0,eps=0)
p5 <- plot_sim2(lambda,os_df_PostDiff_01,c=0.1,eps=0)
p6 <- plot_sim2(lambda,os_df_TTTS_01,c=0,eps=0.1)

lambda <- 0.2
p7 <- plot_sim2(lambda,os_df_TS,c=0,eps=0)
p8 <- plot_sim2(lambda,os_df_PostDiff_01,c=0.1,eps=0)
p9 <- plot_sim2(lambda,os_df_TTTS_01,c=0,eps=0.1)

combined <- p1 + p2+p3+p4+p5+p6+p7+p8+p9 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
#save result into differnt df
#setting 0
c <- 0
eps <- 0
os_df <- os_df_TS

#setting 1
c <- 0.1
eps <- 0
#os_df_PostDiff_01 <- os_df
os_df <- os_df_PostDiff_01

#setting 2
c <- 0
eps <- 0.1
#os_df_TTTS_01 <- os_df
os_df <- os_df_TTTS_01


#
lambda <- 0.1
delta <- c(1:100)/100
score <- (delta/2+lambda*log(2))/(lambda*log(1+exp(delta/lambda)))
plot(delta,score,main = 'lambda=0.1')

lambda <- 0.1
delta <- c(1:100)/100
score <- (delta)/(lambda*log(1+exp(delta/lambda)))
plot(delta,score,main = 'lambda=0.1')



