source("Contextual configration.R")
library(progress)
#note to myself June 27th
#what are some barriers?
#1. Sim/Simulation similar files. what I can copy from sim?
#2. stick to local
#3. how to split the code. One file only for functions, the other only for settings?
#4. there's some fixed settings. How to deal with it?
#5. what else 

# commit message:
#1. need to rename file, revise, and edit AP_prior_sim and Bandit_sim function

AP_prior_sim <- function(B=10000){
  epsilon <- 0 #no need to change, UR ratio
  action_name<<-'a'#should match whatever in the formula defined after this
  beta_name<<-'b'
  true_beta_name<<-'mu'
  context_name<<-'c'
  name_list <<- c(action_name,beta_name,true_beta_name,context_name)
  process_reward_model(true_reward,reward_model,name_list)
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
Bandit_sim <- function(B=10000){
  epsilon <- 0 #no need to change, UR ratio
  action_name <<- 'a'#should match whatever in the formula defined after this
  beta_name <<- 'b'
  true_beta_name<<- 'mu'
  context_name<<- 'c'
  name_list <<-  c(action_name,beta_name,true_beta_name,context_name)
  process_reward_model(true_reward,reward_model,name_list)
  
  bb <- B #number of simulations
  k <- str_count(reward_model,'b') #by default the dim 'k' is the number of 'b_i' in reward model
  bayesian_test <- array(dim=c(bb,k)) #array to save results
  bayesian_est <-  array(dim=c(bb,k)) 
  rewards <- array(dim=c(bb, (1+2*num_act) ))
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
    contextual_matrix <- array(dim=c(n,length(index_list[[context_name]])))
    for (j in 1:length(index_list[[context_name]])){
      contextual_matrix[,j] <- sample(1:length(contextual_distribution_table[,j]),
                                      size = n,
                                      prob = contextual_distribution_table[,j],
                                      replace=T)
      contextual_matrix[,j] <-(contextual_matrix[,j]-1 )/(sum(contextual_distribution_table[,j]>0)-1)
    }
    
    res <- Contextual_Bandit(para_priors,reward_setting,n=n,contextual_matrix=contextual_matrix,epsilon=0,burnin=burnin,batch_size=batch_size)
    process_result(res)
    res1 <- process_result(res)
    rewards[i,] <- as.matrix(res1$reward)
    bayesian_test[i,] <- res1$CI[1,]*res1$CI[3,]>0
    bayesian_est[i,] <- res1$CI[2,]
  }
  reward_output <- data.frame(apply(rewards,2,quantile,probs=c(0.025,0.5,0.975)))
  reward_sd <- apply(rewards,2,sd)
  #colnames(reward_output ) <- c('a1=0','a1=1','all')
  colnames(reward_output) <- colnames(res1$reward)
  output2 <- data.frame(apply(bayesian_test,2,mean))
  output3 <- data.frame(apply(bayesian_est,2,quantile,probs=c(0.025,0.5,0.975)))
  
  return(list(rewards=reward_output,reward_sd=reward_sd,bayesian_test=output2,bayesian_est=output3))
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
#c1: TimeOfDay
true_reward <- 'mu0+mu1*a1+mu2*c1*a1'
reward_model <- 'match' 
mu <- c(0.6,0.1,-0.05) #c(0.3,0.1,-0.05), c(0.6,-0.1,0.05)
contextual_distribution_table <- data.frame(c1=c(0.5,0.5))
mu_prior <- c(0,0,0) #c(0.6,0.1,0), c(0.6,0.2,0)
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




AP_prior_sim(B=100)

res <- Bandit_sim(B=100)
#arm1=0 mean reward
res