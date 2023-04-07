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
    df1 <- rbind(df,expand.grid(algorithm=as.character(i),parameter=parameter_list[[i]]))
  }
  df1$algorithm <- as.character(df1$algorithm )
  df2 <- data.frame(effect=effect_list)#at this point you only have effect size. 
  
  df_settings <- merge(df1,df2)
  df_settings$algo_name <- df_settings$algorithm #name with parameter, so can be used for plots
  df_settings$algo_name[!is.na(df_settings$parameter)] <- paste0(df_settings$algo_name[!is.na(df_settings$parameter)],' (', df_settings$parameter[!is.na(df_settings$parameter)],')')
  return(df_settings) 
}




############## main

#Setting <- 'pa=0.6,pb=0.4,n=785,c=0.1,eps=0'

effect_list <- c(0,0.05,0.1)
horizon_list <- c(100,200,500,785)
lambdas <- c(0.01,0.02,0.04,0.08,0.12,0.16,0.2)
#horizon_max <- max(horizon_list)

parameter_list <- list('TS'=c(NA),
                       'TS PostDiff'=c(0.1,0.2),
                       'TT-TS'=c(0.1,0.5),
                       'UR'=c(NA))


df_settings <- create_df_setting(parameter_list,effect_list)
#process epsilon and c
#can write a better work flow later
df_settings$epsilon <- 0
df_settings$c <- 0
df_settings$epsilon[df_settings$algorithm=='UR'] <- 1
df_settings$epsilon[df_settings$algorithm=='TT-TS'] <- df_settings$parameter[df_settings$algorithm=='TT-TS']

df_settings$c[df_settings$algorithm=='UR'] <- 1
df_settings$c[df_settings$algorithm=='TS PostDiff'] <- df_settings$parameter[df_settings$algorithm=='TS PostDiff']


df_res <- data.frame()
for(i in 1:nrow(df_settings)){
  s <- df_settings[i,]
  para <- list(pa=0.5+s$effect/2,pb=0.5-s$effect/2,n=max(horizon_list),eps=s$epsilon,c=s$c)
  res <- sim2(para,B=10,lambdas=lambdas)
  
}

#use merge


#if length 0, should rplace with NA??
#maybe one column: parameter, instead of one of each
knit_simulation_settings <- function(effect_list,algorithm_list){
  #for simulation compare PostDiff v.s. TT-TS
  
}
ef <- c(0,0.05,0.1,0.2,0.3)
c_list <- c(0,0.1,0.2,1)
eps_list <- c(0,0.1,0.2,0.3)
n_list <- 785

res_df <- expand.grid('effect' = ef,'horizon'=n_list, 'c' = c_list,'epsilon'=eps_list)
res_df$p1 <- 0.5+res_df$effect/2
res_df$p2 <- 0.5-res_df$effect/2
#get rid of cases where both epsilon and c >0.
res_df <- res_df[(res_df$c==0)|(res_df$epsilon==0),]


n <- 1000
c <- 0
eps <- 1
B <- 500
os_df <- array(dim=c(length(ef),n,length(lambdas)))

#a function, loop through p1_v, xxx,n_v,xxxxx. stack all settings
#I guees it's actually easy to go from a stacked result to plot. you just filter it, and 
#specify x-y and group


# a function from several vectors into a dataframe


for (j in 1:length(ef)){
  para <- list(pa=0.5+ef[j]/2,pb=0.5-ef[j]/2,n=n,c=c,eps=eps)
  res <- sim2(para,B=B,lambdas=lambdas)
  os_df[j,,] <- res$objective_scores
}

Setting_summary <- c(Setting_summary,Setting)
para <- eval(parse(text=paste0('list(',Setting,')')))
#sim
#para <- list(pa=0.5,pb=0.5,n=1000,c=0,eps=0)


res <- sim2(para,B=10,lambdas=lambdas)
#for each setting, save result in a spreadsheet (two sub, action and reward). then a summay sheet
sheets <- list("Action" = as.data.frame(t(action_history)), "Reward" = as.data.frame(t(reward_history))) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets, paste0(Setting,".xlsx"),col_names=F)
df <- data.frame(score=c(res$objective_scores),time_step=rep(c(1:para$n),length(lambdas)),lambda=rep(lambdas,each=para$n))
algo <- 'TS PostDiff with c = '
algo_para <- para$c
if(para$eps>0){
  algo_para <- para$eps
  algo='TT-TS with epsilon = '
}

tit <- paste0(algo,algo_para,'; pa=',para$pa,', pb=',para$pb)
ggplot(data=df, aes(x=time_step, y=score, group=lambda, color=lambda)) +
  geom_line()+
  ggtitle(tit)


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



