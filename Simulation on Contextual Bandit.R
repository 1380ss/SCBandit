source("Contextual configration.R")
#note July 20th
# reward calculation
# sampling density: should include e.g. 560, 785 etc
#progress bar
#arm numbering. from 0? UR_sample XXX

#note July 19th
#maybe globale_AP_b or something. A list of global things
#or how to create a df whenever needed? auto figure out dim?
#change sample arm to vector com
#note July 18th
#1. estimation of intercept, main effect etc over time
#reward
#3. allocation probability
#change inverse

#note to myself June 27th
#what are some barriers?
#1. Sim/Simulation similar files. what I can copy from sim?
#2. stick to local
#3. how to split the code. One file only for functions, the other only for settings?
#4. there's some fixed settings. How to deal with it?
#5. what else 

# commit message:
#1. need to rename file, revise, and edit AP_prior_sim and Bandit_sim function

Bandit_sim <- function(B=100,cont_mat_size=10000,AP_track=F,Bayes_track=T,sampling_density=2){
  #simulate in a fixed setting
  k <- str_count(reward_model,'b') #by default the dim 'k' is the number of 'b_i' in reward model
  contextual_matrix_large <- array(dim=c(cont_mat_size,num_context))
  colnames(contextual_matrix_large)<- paste(context_name, index_list[[context_name]], sep="")
  
  magic_for(silent = T,progress = T)
  
  for (j in 1:num_context){
    contextual_matrix_large[,j] <- sample(1:length(contextual_distribution_table[,j]),
                                          size = cont_mat_size,
                                          prob = contextual_distribution_table[,j],
                                          replace=T)
    contextual_matrix_large[,j] <-(contextual_matrix_large[,j]-1 )/(sum(contextual_distribution_table[,j]>0)-1)
  }
  #AP_hist_full <- array(dim=c(n,(1+2*num_context),B))
  #Bayes.est_hist_full <- array(dim=c(n,k,B))
  #coef.est_hist_full <- array(dim=c(n,k,B))
  #p.value_hist_full <- array(dim=c(n,k,B))
  #lm_hist_full <- array(dim=c(n,4,B))

  reward_hist <- c()
  for (i in 1:B){
    #probar$tick()
    set.seed(i)
    contextual_matrix <- contextual_matrix_large[sample(1:cont_mat_size,size=n,replace = T),,drop=F]
    res <- Contextual_Bandit(para_priors,reward_setting,n=n,contextual_matrix=contextual_matrix,burnin=burnin,batch_size=batch_size)
    if (AP_track){
      AP_hist_full[,,i] <- get_AP_hist(prior=para_priors,contextual_matrix_large,n=n,X=res$X,y=res$reward_hist,AP_b=1000,sampling_density=sampling_density)
    }
    if (Bayes_track){
      Bayes.est_hist_full[,,i] <- get_Bayes_hist(res$X,res$reward_hist,para_priors,sampling_density=sampling_density)
    }
    #re <- get_lm_hist(res$X,res$reward_hist)
    #coef.est_hist_full[,,i] <- re$est.coef_hist
    #p.value_hist_full[,,i] <- re$p.value_hist
    #reward_hist[i] <- mean(res$reward_hist)
    put(mean(res$reward_hist))
  }
  #return(list(Bayes.est_hist=Bayes.est_hist_full,AP_hist=AP_hist_full,coef.est_hist=coef.est_hist_full,p.value_hist=p.value_hist_full,reward_hist=reward_hist))
  return(magic_result())
}

stepwise_AP <- function(para_post,contextual_matrix_large,AP_b,context_var='c2'){
  #this only works for normal Bayesian model
  #maybe use Stan in the future
  #how to make contextual matrix more standardized?
  #change input to hist! now there's too many, and repeated input
  
  cont_temp <- contextual_matrix_large[sample(1:cont_mat_size,size=AP_b,replace = T),,drop=F]
  act_temp <- as.matrix(array(dim=c(AP_b,num_act)))
  colnames(act_temp) <- paste(action_name, index_list[[action_name]], sep="")
  for (j in 1:AP_b){
    act_temp[j,] <- sample_arms(para_post,num_act,context=cont_temp[j,])
  }
  df <- as.data.frame(cbind(act_temp,cont_temp))
  
  #for i in 1:numcontext
  
  return(c(summarize(group_by_at(df,context_var),a1=mean(a1))[,2][[1]],mean(act_temp)))
  #AP_hist[i,] <- c(summarize(group_by(df,c1),mean=mean(a1))[,2][[1]],mean(act_temp)) #apply(act_temp,2,mean)
}

stepwise_lm <- function(X,y,i){
  det_X <- determinant(t(X[1:i,])%*%(X[1:i,]))
  coef <- NA
  pvalue <- NA
  if((det_X$modulus[1] > -3)&(i>1)){
    lm_setpwise <- summary(lm(y[1:i]~X[1:i,-1]))
    coef <- lm_setpwise$coefficients[,1]
    pvalue <- lm_setpwise$coefficients[,4]
  }
  return(list(coef=coef,pvalue=pvalue))
}

get_Bayes_hist <- function(X,y,para_priors,sampling_density=2){
  #sampling_density is roughly 'how often you do a sample', it's a log rate anyway (sample less and less)
  n <- dim(X)[1]
  Bayes.est_hist <- array(dim=c(n,dim(X)[2]))
  i <- 1
  imp <- 0 #importance, accumulate as time increases
  while(i<=n){
    imp <- imp+1/i
    if( imp>1/sampling_density){
      re <- blinear_update(X[1:i,,drop=F],y[1:i,drop=F],para_priors)
      Bayes.est_hist[i,] <- re$mu
      imp <- imp-1/sampling_density
    }
    i <- i+1
  }
  return(Bayes.est_hist)
}
get_lm_hist <- function(X,y,sampling_density=2){
  n <- dim(X)[1]
  est.coef_hist <- array(dim=c(n,dim(X)[2]))
  p.value_hist <- array(dim=c(n,dim(X)[2]))
  i <- 1
  imp <- 0 #importance, accumulate as time increases
  while(i<=n){
    imp <- imp+1/i
    if( imp>1/sampling_density){
      re <- stepwise_lm(X,y,i)
      est.coef_hist[i,] <- re$coef
      p.value_hist[i,] <- re$pvalue
      imp <- imp- 1/sampling_density
    }
    
    i <- i+1

  }
  
  return(list(est.coef_hist=est.coef_hist,p.value_hist=p.value_hist))
}
get_AP_hist <- function(priors,contextual_matrix_large,n,X,y,AP_b=3000,sampling_density=2){
  #this only works for normal Bayesian model
  #maybe use Stan in the future
  #how to make contextual matrix more standardized?
  
  #force match dim
  para_post <- blinear_update(X[1,,drop=F],y[1,drop=F],priors=priors)
  AP_hist <- array(dim=c(n,length(stepwise_AP(para_post,contextual_matrix_large,AP_b))))
  i <- 1
  imp <- 0 #importance, accumulate as time increases
  while(i<=n){
    
    imp <- imp+1/i
    if( imp>1/sampling_density){
      para_post <- blinear_update(X[1:i,,drop=F],y[1:i,drop=F],priors=priors)
      AP_hist[i,] <- stepwise_AP(para_post,contextual_matrix_large,AP_b)
      imp <- imp- 1/sampling_density
    }
    
    i <- i+1
    
  }
  
  return(AP_hist)
}



Bandit_sim_unisetting <- function(bb=10000){
  #simulate in a fixed setting
  k <- str_count(reward_model,'b') #by default the dim 'k' is the number of 'b_i' in reward model
  bayesian_test <- array(dim=c(bb,k)) #array to save results
  bayesian_est <-  array(dim=c(bb,k)) 
  rewards <- array(dim=c(bb, (1+2*num_context) ))
  ap_df <- array(dim=c(bb,n,3 ))
  probar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                             total = bb,
                             complete = "=",   # Completion bar character
                             incomplete = "-", # Incomplete bar character
                             current = ">",    # Current bar character
                             clear = FALSE,    # If TRUE, clears the bar when finish
                             width = 100)      # Width of the progress bar
  
  for (i in 1:bb){
    probar$tick()
    set.seed(i)
    #n sample size
    #k number of betas
    
    #need to change contextual matrix each time
    res <- Contextual_Bandit(para_priors,reward_setting,n=n,contextual_matrix=contextual_matrix,epsilon=0,burnin=burnin,batch_size=batch_size,record_AP=T)
    
    
    
    process_result(res)
    res1 <- process_result(res)
    rewards[i,] <- as.matrix(res1$reward)
    bayesian_test[i,] <- res1$CI[1,]*res1$CI[3,]>0
    bayesian_est[i,] <- res1$CI[2,]
    ap_df[i,,] <- res$AP_hist
  }
  reward_output <- data.frame(apply(rewards,2,quantile,probs=c(0.025,0.5,0.975)))
  reward_sd <- apply(rewards,2,sd)
  #colnames(reward_output) <- c('a1=0','a1=1','all')
  colnames(reward_output) <- colnames(res1$reward)
  output2 <- data.frame(apply(bayesian_test,2,mean))
  output3 <- data.frame(apply(bayesian_est,2,quantile,probs=c(0.025,0.5,0.975)))
  
  return(list(rewards=reward_output,reward_sd=reward_sd,bayesian_test=output2,bayesian_est=output3,ap_df=ap_df))
}



AP_prior_sim <- function(B=10000,epsilon=0){
  para_post <- para_priors
  
  contextual_matrix <- array(dim=c(B,num_context))
  act_hist <- as.data.frame(array(NA,dim=c(B,num_act)))
  for (j in 1:num_context){
    contextual_matrix[,j] <- sample(1:length(contextual_distribution_table[,j]),
                                    size = B,
                                    prob = contextual_distribution_table[,j],
                                    replace=T)
    contextual_matrix[,j] <-(contextual_matrix[,j]-1 )/(sum(contextual_distribution_table[,j]>0)-1)
  }
  for (i in 1:B){
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
  }
  if(dim(act_hist)[2]>1){
    return(apply(act_hist,2,mean))#allocation to arm1 (GetRationale)
  } else{
    return(mean(act_hist[,1]))#allocation to arm1 (GetRationale)
  }
  
}

##############   list of tunning parameters
#1. 'n', sample size,
#2. 'contextual_matrix_large', a matrix describing the dsitribution of contexts. Each columns denotes a differnt context (by default it's c1, c2,...). Can modify to be more flexible
#3. 'mu': true beta parameters generating rewards;
#4. name_list <- c(action_name,beta_name,context_name). Name identifiers
#5. true_reward <- 'b0+b1*a1+b2*c1+b3*a1*c1'
#   reward_model <- 'b0+b1*a1+b2*c1+b3*a1*c1'
#6. 'sig_err', variance of reward error
#7. 'priors', (a0=2,b0=1,mu0=rep(0,k),L0=0.01*diag(k))

n <- 560
batch_size <- 4
burnin <- 8

############# study B2 
##################

#Rating = intercept + BetaRationale*GetRationale + BetaRationaleTimeOfDay*GetRationale*TimeOfDay
#b0: Intercept
#b1: BetaRationale
#b2: BetaRationaleTimeOfDay
#a1: GetRationale
#a2: NoRationale
#c1: TimeOfDay

#DidActivity = intercept + BetaRationale*GetRationale + BetaRationaleTimeOfDay*GetRationale*TimeOfDay
#b0: Intercept
#b1: BetaRationale
#b2: BetaRationaleTimeOfDay
#a1: GetRationale
#c1: TimeOfDay
true_reward <- 'mu0+mu1*a1+mu2*c1+mu3*c2*a1'
reward_model <- 'match' 
mu <- c(0.4,0.1,0.1,-0.15) #c(0.3,0.1,-0.05), c(0.6,-0.1,0.05)
contextual_distribution_table <- data.frame(c1=c(0.1,0.2,0.3,0.4),
                                            c2=c(0.4,0.1,0.1,0.4))
mu_prior <- c(0,0,0,0) #c(0.6,0.1,0), c(0.6,0.2,0)
var_diag <- 5 #inverse.  500, 0.05


#study C1
#true_reward <- 'mu0+mu1*a1+mu2*a2+mu3*c1+mu4*a1*c1+mu5*a2*c1'

#Study C2
#true_reward <- 'mu0+mu1*a1+mu2*a2+mu3*a3+mu4*a1*c1+mu5*a2*c1+mu6*a3*c1'

#Study C3
#true_reward <- 'mu0+mu1*a1+mu2*c1*a1'


#test: can we get rid of it??
sig_err <- 1/36 #variance of reward error
reward_setting <- list(mu=mu,sig_err=sig_err)
reward_setting$reward_generation <- function(mean_reward){
  #y <- median(c(0,1,round(y*4)/4))
  #y <- rnorm(1,mean=mean_reward,sd=sqrt(sig_err))
  #y <- median(c(0,1,y))
  y <- rbernoulli(1,p=mean_reward)
  return(y)
}


############ fixed setting
#######################

if (reward_model=='match'){
  reward_model <- paste(str_split(true_reward,'mu')[[1]],collapse = 'b')
}

num_context <- length(unique(get_indices(paste0(true_reward,reward_model),'c')))



Var_cov_prior <- var_diag*diag(length(mu_prior))
para_priors <- list(a=2,b=1,mu=mu_prior,L=Var_cov_prior)
process_reward_model(true_reward,reward_model,name_list)



res <- Bandit_sim(B=100,AP_track = F, Bayes_track = F,sampling_density=2)








#res_var5 <- res
res <- res_var5$AP_hist
i <- 5
df <- apply(res[,i,],1,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
ind <- !is.na(df[1,])
df <- df[,ind]
plot((1:n)[ind],df[2,],type='l',
     main='Average allocation probability (across all) \n over time with 95% CI',
     xlab='Time Step',
     ylab='Allocation Prob',
     ylim=c(0,1))
lines((1:n)[ind],df[1,],type='l',col='red')
lines((1:n)[ind],df[3,],type='l',col='red')























instance=1
i <- 2
x <- !is.na(res$Bayes.est_hist[,1,1])
plot(c(1:n)[x],res$Bayes.est_hist[x,i,instance],type='l',
     ylim=c(0,0.5),
     xlab = 'time step',
     ylab='estimation of main effect (Bayes)',
     main = paste0('Bayes Estimation of main effect over time, \nwith true main effect =',mu[2], ', variance =',1/var_diag))
abline(h=0.2,col='red')

















#res_ef01 <- res
res_ef005 <- res
ind <- !is.na(re01)

res <- res_ef01
re01 <- apply(res$p.value_hist<0.05,c(1,2),mean)
res <- res_ef005
re005 <- apply(res$p.value_hist<0.05,c(1,2),mean)

write.csv(cbind(c(1:n)[ind],re01[ind,2],re01[ind,4],re005[ind,2],re005[ind,4]),file='~/re.csv')


res_var005 <- res
#res_var500 <- res
#res_var5 <- res
#then, change var to 0.5, and 500
instance <- 2

res <- res_var005
i <- 2
x <- !is.na(res$lm_hist[,1,1])
plot(c(1:n)[x],res$lm_hist[x,i,instance],type='l',
     xlab = 'time step',
     ylab='estimation of main effect',
     main = 'Estimation of main effect over time, \nwith true main effect =0.1, variance = 20')
abline(h=0.1,col='red')

res <- res_var5
i <- 2
x <- !is.na(res$lm_hist[,1,1])
plot(c(1:n)[x],res$lm_hist[x,i,instance],type='l',
     xlab = 'time step',
     ylab='estimation of main effect',
     main = 'Estimation of main effect over time, \nwith true main effect =0.1, variance = 0.2')
abline(h=0.1,col='red')

res <- res_var500
i <- 2
x <- !is.na(res$lm_hist[,1,1])
plot(c(1:n)[x],res$lm_hist[x,i,instance],type='l',
     xlab = 'time step',
     ylab='estimation of main effect',
     main = 'Estimation of main effect over time, \nwith true main effect =0.1, variance = 0.002')
abline(h=0.1,col='red')









#AP_b <- 1000
contextual_matrix <- as.matrix(array(dim=c(AP_b,length(index_list[[context_name]]))))
heading <- c()
for (b in 1:length(index_list[[context_name]])){
  heading <- paste0('c',b)
}
colnames(contextual_matrix) <- heading
for (j in 1:length(index_list[[context_name]])){
  contextual_matrix[,j] <- sample(1:length(contextual_distribution_table[,j]),
                                  size = AP_b,
                                  prob = contextual_distribution_table[,j],
                                  replace=T)
  contextual_matrix[,j] <-(contextual_matrix[,j]-1 )/(sum(contextual_distribution_table[,j]>0)-1)
}























######################## old
######
res_var5 <- res
res <- res_var5$AP_hist
i <- 5
df <- apply(res[,i,],1,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
ind <- !is.na(df[1,])
df <- df[,ind]
plot((1:n)[ind],df[2,],type='l',
     main='Average allocation probability (all group) over time with 95% CI',
     xlab='Time Step',
     ylab='Allocation Prob',
     ylim=c(0,1))
lines((1:n)[ind],df[1,],type='l',col='red')
lines((1:n)[ind],df[3,],type='l',col='red')



i <- 1
df <- apply(res[,i,],1,quantile,probs=c(0.025,0.5,0.975))
plot(1:n,df[2,],type='l',
     main='Average allocation probability (c1=0) over time with 95% CI',
     xlab='Time Step',
     ylab='Allocation Prob',
     ylim=c(0,1))
lines(1:n,df[1,],type='l',col='red')
lines(1:n,df[3,],type='l',col='red')


i <- 2
df <- apply(res[,i,],1,quantile,probs=c(0.025,0.5,0.975))
plot(1:n,df[2,],type='l',
     main='Average allocation probability (c1=1) over time with 95% CI',
     xlab='Time Step',
     ylab='Allocation Prob',
     ylim=c(0,1))
lines(1:n,df[1,],type='l',col='red')
lines(1:n,df[3,],type='l',col='red')


res1 <- Contextual_Bandit(para_priors,reward_setting,n=100,contextual_matrix,burnin,batch_size,epsilon,record_AP=T,AP_b=AP_b,record_coef=T)

#seom code create that csv
#cahnge matrix colnames names
#X, action, contextual, reward,
output_csv <- cbind(res1$X,res1$his,res1$AP_hist,res1$lm_coef_hist)
write.csv(output_csv,file='~/MHA_stepwise_output.csv')

df1 <- array(dim=c(2000,3))
for (i in 1:20){
  df1[((i-1)*100+1):(i*100),1] <- res$ap_df[i,,3]
  df1[((i-1)*100+1):(i*100),2] <- i
  df1[((i-1)*100+1):(i*100),3] <- 1:100
}

df1 <- as.data.frame(df1)
colnames(df1) <- c('Allocation_probability','trial','Time_Step')
ggplot(data=df1,aes(x=Time_Step,y=Allocation_probability,group=trial))+
  geom_line()
  ggtitle('allocation prob over time')


df <- cbind(c(1:n),res$AP_hist,res$lm_coef_hist)
df <- df[!is.na(df[,4]),]

#ap
plot(df[,1],df[,2],main = 'Allocation Probability to arm1',
     xlab = 'time step',
     ylab = 'Allocation probability')

plot(df[,1],df[,4],main = 'estimation of intercetp over time (LM)',
     xlab = 'time step',
     ylab = 'estimation of intercept')

AP_prior_sim(B=100)

res <- Bandit_sim_unisetting(B=20)
#arm1=0 mean reward
res
