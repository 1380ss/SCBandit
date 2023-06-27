library(invgamma)
library(pracma)
library(stringr)
library(mvtnorm)
library(tidyverse)
library(readr)
library(rebus)
library(dplyr)
options(dplyr.summarise.inform = FALSE)
#read some scripts

#standardize your inpout/output
#also, if output is standardized, you can use it without too much review!


#We used this R code to simulate:
# 1. in three arm case, performance of Uniform Random v.s. Top two TS
# 2. in two arm, one contextual varaible (high/low accuracy) case, perfroamce of TS-contextual v.s. TT-TS v.s. UR v.s. UR-then-exploit



################### minor functions ###########################################
get_indices <- function(expr,var_name,uniq=F){
  #get index numbers for a variable in a string
  #expr is the full string of the expression (e.g. reward model etc)
  ## example: get_indices('b0+b1*a1+b2*c1+b3*a1*c1','c')
  indices <- parse_number(str_split(expr,char_class(var_name))[[1]])[-1]
  indices <- indices[!is.na(indices)]
  if(uniq){
    #get unique, ordered indices
    indices <- unique(indices)
    indices <- indices[order(indices)]
  }
  return(indices)
}

get_IndexMatrix <- function(expr,index_list,sep='+'){
  ## example:get_IndexMatrix('b0+b1*a1+b2*c1+b3*a1*c1')
  s <- str_split(expr,char_class(sep))[[1]]
  m <- matrix(nrow=length(index_list),ncol=length(s))
  row.names(m) <-names(index_list) 
  for (i in name_list){
    for (j in 1:length(s)){
      ind <- parse_number(str_split(s[j],i)[[1]])[-1]
      if(length(ind>0)){
        m[i,j] <- ind
      }
    }
  }
  return(m)
}

index_mapping <- function(im,index_list){
  for(i in row.names(im)){
    for (j in 1:ncol(im)){
      if(!is.na(im[i,j])){
        im[i,j] <- which(index_list[[i]]==im[i,j])
      }
    }
  }
  return(im)
}

m2v <- function(m,v){
  #m is the matrix
  #v is the value list
  #will skip NAN values.
  #can also be viewed as the reverse of 'index_mapping' function
  for(i in row.names(m)){
    for (j in 1:ncol(m)){
      if(!is.na(m[i,j])){
        m[i,j] <- v[[i]][m[i,j]]
      }
    }
  }
  return(m)
}


#################### end of minor functions #################################

UR_sample_actions <- function(num_act,onehot=T){
  if (onehot){
    action_draw <- sample(1:(1+num_act),size=1)
    action_vector <- rep(0,num_act) #need to confirm same length as actual action array
    if (action_draw<=num_act){
      action_vector[action_draw] <- 1
    }
  } else {
    stop("need to re-work on UR_sample_actions function, no onehot=F option")
  }
  return(action_vector=action_vector)
}

blinear_update <- function(X,y,priors){
  a0 <- priors$a
  b0 <- priors$b
  mu0 <- priors$mu
  L0 <- priors$L
  n <- dim(X)[1]
  H <- t(X)%*%X
  #b_hat <- inv(H)%*%(t(X)%*%y)
  Ln <- H+L0
  mun <- inv(Ln) %*% ((t(X)%*%y)+L0%*%mu0)
  an <- a0+n/2
  bn <- b0+(t(y)%*%y+t(mu0)%*%L0%*%mu0-t(mun)%*%Ln%*%mun)/2
  return(list(a=an,b=bn,mu=mun,L=Ln))
}

#I think each action variable should be binary. Can explain more.
mat2exp <- function(mat){
  #input: a matrix, with rows corresponding to variable names, and entries corresponding to indices
  #output: a text of expression of everything collasped
  #return 1 if everything is NA.
  var_names <- row.names(mat)
  formula_pieces <- c()
  mat_copy <- mat #I think that may lead to bugs.. Better if just create an NA matrix with same row names
  for(i in var_names){
    for(j in 1:ncol(mat)){
      if(!is.na(mat[i,j])){
        mat_copy[i,j] <- paste0(i,'[',mat[i,j],']')#can do -1 for Python
      }
    }
  }
  for(j in 1:ncol(mat)){
    formula_pieces[j] <-paste(mat_copy[!is.na(mat_copy[,j]),j],collapse = '*') 
  }
  full_formula <- paste(formula_pieces,collapse='+')
  if(str_count(full_formula)==0){
    full_formula <- '1'
  }
  return(full_formula)
}

mat2coef <- function(mat,para_name){
  #mat2coef(regression_im_sort_ind,'a')
  #input: 'mat': a matrix, with rows corresponding to variable names, and entries corresponding to indices
  #input: 'para_name': the parameter you look for coef
  #output: a vector of texts of coefs
  coefs <- c()
  for (i in unique(mat[para_name,])){
    if(!is.na(i)){
      coefs[i] <- mat2exp(mat[rownames(mat) !=para_name,which(mat[para_name,]==i),drop = F])
    }
  }
  return(coefs)
}

process_reward_model <- function(true_reward,reward_model,name_list){
  #process_reward_model(true_reward,reward_model,name_list)
  #produce action_indicator function, input must be b= and c=
  #produce reward_eval function, input must be a, c and mu
  #produce get_design_matrix 
  index_list <<- list()
  for (i in name_list){
    index_list[[i]] <<- get_indices(paste0(true_reward,reward_model),i,uniq=T)
  }
  
  #act_indices <- get_indices(paste0(true_reward,reward_model),action_name,uniq=T)
  num_act <<- length(index_list[[action_name]])
  
  #beta_indices <- get_indices(paste0(true_reward,reward_model),beta_name,uniq=T)
  num_beta <<- length(index_list[[beta_name]])
  
  regression_im_ori_ind <- get_IndexMatrix(reward_model,index_list) #align with original index in formula
  regression_im_sort_ind <- index_mapping(regression_im_ori_ind,index_list)
  
  true_im_ori_ind <- get_IndexMatrix(true_reward,index_list) #align with original index in formula
  true_im_sort_ind <- index_mapping(true_im_ori_ind,index_list)
  
  coefs <- mat2coef(regression_im_sort_ind,action_name)
  reward_eval_text <- mat2exp(true_im_sort_ind)
  dm <- mat2coef(regression_im_sort_ind,beta_name)
  eval(parse(text= paste0('action_indicator<<-function(b,c,onehot=T){ ifelse(onehot, order(-c(',paste(coefs,collapse=','),',0))[1] ,  c(',paste(coefs,collapse=','),')>0) }') ))
  eval(parse(text= paste0('reward_eval<<-function(a,c,mu){return(',reward_eval_text,')}') ))
  eval(parse(text= paste0('get_design_matrix<<-function(a,c){c(',paste(dm,collapse=','),')>0}') ))
  return(list(action_indicator=coefs,get_design_matrix=dm,reward_eval=reward_eval_text))
  #return(list(num_act=num_act,coef_list=coef_list,act_list=act_list))
}

Contextual_Bandit <- function(para_priors,reward_setting,n,contextual_matrix,burnin,batch_size,epsilon,record_AP=F){
  #it's quite clear what this does if you look at return. 
  if (dim(contextual_matrix)[2]!=length(index_list[[context_name]])){
    stop("Incorrect number of columns in contextual_matrix")
  }
  
  para_post <- para_priors
  
  act_hist <- array(NA,dim=c(n,num_act))
  #allocation prababilty 
  
  #action 00 -> arm 1
  #action 01 ->arm2
  
  
  #AP_hist <- array(NA, dim=c(n,(num_act)+1))
  #column 1-> arm1 ect...
  
  X <- matrix(nrow=n,ncol=length(para_priors$mu))
  y <- matrix(nrow=n,ncol=1)
  value_list <- list()
  
  #burinin phase
  if(burnin>0){
    for (i in 1:burnin){
      #notes:
      #for arm, if we have k arms, than we need k-1 action variable. (think of regression equation).
      #their sum need to be <=1
      
      #original UR
      #act_hist[i,] <- 1*(runif(num_act)>1/2)
      
      #new design, sum <=1
      act_hist[i,] <- UR_sample_actions(num_act = num_act)
      
      
      
      #three arm_> 2 action variable, 11 resample 
      X[i,] <- 1*get_design_matrix(a=act_hist[i,],c=contextual_matrix[i,])
      mean_reward <- reward_eval(a=act_hist[i,],c=contextual_matrix[i,],mu=reward_setting$mu)
      y[i] <- reward_setting$reward_generation(mean_reward)
      
    }
  }
  
  #should I do update here?
  para_post <- blinear_update(X[1:i,,drop = F],y[1:i],para_priors)
  batch_ind <- 0
  if(burnin<n){
    for (i in (burnin+1):n){
      batch_ind <- batch_ind+1
      #sample beta and sigma
      sigma_n <- rinvgamma(1,para_post$a,para_post$b) #sigma here is actually sigma square
      beta_n <- rmvnorm(1,mean=para_post$mu,sigma = sigma_n*inv(para_post$L))
      
      if(runif(1)<epsilon){
        act_hist[i,] <- UR_sample_actions(num_act = num_act)
      } else {
        #if onehot =T, only return a single arm each time. Otherwise, a vector...
        #someone can change it
        arm_ind <- action_indicator(b=beta_n,c=contextual_matrix[i,],onehot=T)
        
        act_hist[i,] <- 0
        if (arm_ind<=num_act){
          act_hist[i,arm_ind] <- 1
        }
        
        
        
      }
      
      X[i,] <- 1*get_design_matrix(a=act_hist[i,],c=contextual_matrix[i,])
      mean_reward <- reward_eval(a=act_hist[i,],c=contextual_matrix[i,],mu=reward_setting$mu)
      y[i] <- reward_setting$reward_generation(mean_reward)
      #
      if (batch_ind==batch_size){
        batch_ind <- 0
        para_post <- blinear_update(X[1:i,,drop = F],y[1:i],para_priors)
      }
      
      
    }
  }
  
  para_post <- blinear_update(X[1:i,,drop = F],y[1:i],para_priors)
  
  colnames(act_hist) <- paste(action_name, index_list[[action_name]], sep="")
  colnames(contextual_matrix)<- paste(context_name, index_list[[context_name]], sep="")
  colnames(y) <- 'reward'
  his <- cbind(act_hist,contextual_matrix,y)
  
  return(list(para_post=para_post,X=X,y=y,his=his))
}

process_result <- function(res,B=1000,inds,names){
  y <- res$y
  betas <- array(dim=c(B,length(res$para_post$mu)))
  his <- as_tibble(res$his)
  col <- c('all')
  reward <- mean(y)
  for(i in 1:num_context){
    col <- c(col,paste0('c',i,'=0'),paste0('c',i,'=1'))
    reward <- c(reward ,summarize(group_by_at(his,paste0('c',i)),reward=mean(reward))$reward)
  }
  reward <- data.frame(t(reward))
  colnames(reward) <- col
  mse <- summarize(group_by(his,a1),mse=sd(reward)/sqrt(n()))$mse
  for (i in 1:B){
    sigma_n <- rinvgamma(1,res$para_post$a,res$para_post$b) #sigma here is actually sigma square
    betas[i,] <- rmvnorm(1,mean=res$para_post$mu,sigma = sigma_n*inv(res$para_post$L))
  }
  
  return( list(reward=reward,
               MSE=c(mse,sd(y)/sqrt(length(y))),
               CI=apply(betas,2,quantile,probs=c(0.025,0.5,0.975))))
  
  #betas <- array(dim=c(B,length(mun)))
  #his <- as_tibble(his)
  #his_sum <- summarize(group_by(his,a1,c1),p_hat=mean(reward),sd=sd(reward),n=n())
  #Wald <- Wald_for_cont(his_sum)
  #his_sum <- summarize(group_by(his,a1),reward=mean(reward))
  #for (i in 1:B){
  #  sigma_n <- rinvgamma(1,an,bn) #sigma here is actually sigma square
  #  betas[i,] <- rmvnorm(1,mean=mun,sigma = sigma_n*inv(Ln))
  #}
  #nx <- dim(X)[2]-1
  #ind <- 2^(c(1:nx)-1)
  #xy <-as_tibble(cbind(X,y))
  #colnames(xy)[inds] <- names
  #X_group <- summarise(group_by_at(xy,inds),m=mean(y),con=sum(y)/mean(y))
  #return(list(CI=apply(betas,2,quantile,probs=c(0.025,0.5,0.975)),
  #            reward=X_group,Wald=Wald))
}
