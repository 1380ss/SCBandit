source("Contextual Bandit Algorithm.R")


######################################################################
##############   everything that needs to be defined, begin  ############################
#1. 'n', sample size,
#2. 'contextual_matrix_large', a matrix describing the dsitribution of contexts. Each columns denotes a differnt context (by default it's c1, c2,...). Can modify to be more flexible
#3. 'mu': true beta parameters generating rewards;
#4. name_list <- c(action_name,beta_name,context_name). Name identifiers
#5. true_reward <- 'b0+b1*a1+b2*c1+b3*a1*c1'
#   reward_model <- 'b0+b1*a1+b2*c1+b3*a1*c1'
#6. 'sig_err', variance of reward error
#7. 'priors', (a0=2,b0=1,mu0=rep(0,k),L0=0.01*diag(k))

n <- 1000 # 100 #number of participants
B <- 10000 #used below. B needs to be large enough so you get a good approximation of the distribution of contexts
nc <- 2 #total number of context involved
prob_high <- 0.5
#need a distribution of context
#here I just make them randomly normal. But you can generate any distribution based on your assumption.
#It's not necessary to make it a large matrix. e.g. if you exhaustively list all possible combination of contexts and reflect their true distribution probability, it's even better.
contextual_matrix_large <- array(rnorm(B*nc),dim=c(B,nc)) #generate the context for all participants
contextual_matrix_large[,1] <- rbinom(B,1,prob=prob_high)

#define true reward equation using a text string in the following format
#b_i is beta parameter
#c_i is contextual varaibles (e.g. )
#a_i is actions (e.g. link v.s. no-link)

#I think currectly need to fix using a,b and c. 
#numbers better start at 1, and increase by 1.
action_name='a'#should match whatever in the formula defined after this
beta_name='b'
true_beta_name='mu'
context_name='c'
name_list <- c(action_name,beta_name,true_beta_name,context_name)
true_reward <- 'mu0+mu1*a1+mu2*c1+mu3*a1*c1'
reward_model <- 'b0+b1*a1+b2*c1+b3*a1*c1'
reward_model <- 'b0+b1*a1+b2*a1*c1'
sig_err <- 1/36 #variance of reward error
#true_beta <- c(b0,b1,b2,b3) #true value of b_i in true_reward 
mu <- c(0.5,3/8,-1/4,-5/8) #true value of mu_i in true_reward 
reward_setting <- list(mu=mu,sig_err=sig_err)

reward_setting$reward_generation <- function(mean_reward){
  y <- rnorm(1,mean=mean_reward,sd=sqrt(sig_err))
  y <- median(c(0,1,round(y*4)/4))
  return(y)
}
#true_reward <- 'b0+b1*a1'
#reward_model <- 'b0+b1*a1'
#sig_err <- 1/36 #variance of reward error
#true_beta <- c(0.5,0.125) #true value of b_i in true_reward 
#don't need to define k below, but should choose priors
k <- str_count(reward_model,'b') #by default the dim 'k' is the number of 'b_i' in reward model
para_priors <- list(a=2,b=1,mu=rep(0,k),L=0.01*diag(k))

#checklist of all you need. Not sure if we use this actually...
sim_setting <- list(n=n,
                    name_list=name_list,
                    true_beta=true_beta,
                    true_reward=true_reward,
                    reward_model=reward_model,
                    sig_err=sig_err,
                    priors=priors,
                    contextual_matrix_large=contextual_matrix_large)

#########################   end  ############################



process_reward_model(true_reward,reward_model,name_list)


bb <- 300
bt <- array(dim=c(bb,k))
rewards <- array(dim=c(bb,3))
cont <- array(dim=c(bb,k))
w <- array(dim=c(bb,3))

probar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = bb)      # Width of the progress bar
for (i in 1:bb){
  probar$tick()
  set.seed(i)
  #n sample size
  #k number of betas
  res <- Contextual_Bandit(para_priors,reward_setting,n=100,contextual_matrix,epsilon=0,burnin=100,batch_size=1)
  res1 <- process_result(res)
  rewards[i,] <- res1$reward
  #res0
  bt[i,] <- res1$CI[1,]*res1$CI[3,]>0
  #rew[i,] <- pull(res0$reward[,3])
  #cont[i,] <- pull(res0$reward[,4])
  # w[i,] <- res0$Wald
  #summary(lm(res$y~res$X))
  #summary(glm(res$y~res$X,family=binomial))
}

apply(rewards,2,mean)
apply(bt,2,mean)