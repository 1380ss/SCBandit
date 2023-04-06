#Goal: givne output, 
#load that function TS2b
#1. function calculating optimization of objective with differnt lambda
# look at your overleaf..
#2. calculate empirical performance (at each step, each 'c' and epsilon.)
library(progress)
#https://r-coder.com/progress-bar-r/ 

library(ggplot2)
library(patchwork)
source("Algorithm for two arm.R")
source("Algorithms for general use.R")

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
entropy_objective_optimal_2arm(lambda=0.2,delta=0.2,x=0.9)

#sim
para <- list(pa=0.5,pb=0.5,n=1000,c=0,eps=0)

lambdas <- c(0.05,0.1,0.2,0.3,0.5,1,2,3)
res <- sim2(para,B=1000,lambdas=lambdas)
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



