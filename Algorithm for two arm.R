#updated April 3rd 2024

#I shouldn't have library loaded here...

#TS algorithm for 2-arm, binary reward case. Can use PostDiff and epsilon-UR
TS2b <- function(pa,pb,n,c,eps,burnin=0,ap_output=F){
  # n is the sample size (785 or 197 in our case)
  arm_a_successes <- c(0)
  arm_a_failures <- c(0)
  arm_b_successes <- c(0)
  arm_b_failures <- c(0)
  action_each_step <- c()
  reward_each_step <- c()
  arm_a_s <- c(0)
  arm_a_f <- c(0)
  arm_b_s <- c(0)
  arm_b_f <- c(0)
  i <- 0
  
  while (i < burnin){
    i <- i+1
    draw_a <- runif(1)
    draw_b <- runif(1)
    action_each_step[i] <- (draw_a>draw_b)
    if (action_each_step[i]){
      reward_each_step[i] <- (runif(1)<pa)
      arm_a_s <- arm_a_s+reward_each_step[i]
      arm_a_f <- arm_a_f+1-reward_each_step[i]
    } else{
      reward_each_step[i] <- (runif(1)<pb)
      arm_b_s <- arm_b_s+reward_each_step[i]
      arm_b_f <- arm_b_f+1-reward_each_step[i]
    }
    arm_a_successes[i+1] <- arm_a_s
    arm_a_failures[i+1] <- arm_a_f
    arm_b_successes[i+1] <- arm_b_s
    arm_b_failures[i+1] <- arm_b_f
  }
  
  while (i < n){
    i <- i+1
    if (runif(1)<eps){
      #epsilon-UR step
      draw_a <- runif(1)
      draw_b <- runif(1)
    } else {
      #TS-step
      draw_a <- rbeta(1,arm_a_successes[i]+1,arm_a_failures[i]+1)
      draw_b <- rbeta(1,arm_b_successes[i]+1,arm_b_failures[i]+1)
      if (abs(draw_a-draw_b)<c){
        #PostDiff step
        draw_a <- runif(1)
        draw_b <- runif(1)
      } else {
        draw_a <- rbeta(1,arm_a_successes[i]+1,arm_a_failures[i]+1)
        draw_b <- rbeta(1,arm_b_successes[i]+1,arm_b_failures[i]+1)
      }
    }
    action_each_step[i] <- (draw_a>draw_b)
    if (action_each_step[i]){
      reward_each_step[i] <- (runif(1)<pa)
      arm_a_s <- arm_a_s+reward_each_step[i]
      arm_a_f <- arm_a_f+1-reward_each_step[i]
    } else{
      reward_each_step[i] <- (runif(1)<pb)
      arm_b_s <- arm_b_s+reward_each_step[i]
      arm_b_f <- arm_b_f+1-reward_each_step[i]
    }
    arm_a_successes[i+1] <- arm_a_s
    arm_a_failures[i+1] <- arm_a_f
    arm_b_successes[i+1] <- arm_b_s
    arm_b_failures[i+1] <- arm_b_f
  }
  #na: number of times arm A was pulled. nb: same
  
  arm_a_successes <- arm_a_successes[-1]
  arm_a_failures <- arm_a_failures[-1]
  arm_b_successes <- arm_b_successes[-1]
  arm_b_failures <- arm_b_failures[-1]
  
  na <- arm_a_successes+arm_a_failures
  nb <- arm_b_successes+arm_b_failures
  pa_est <- arm_a_successes/na
  pb_est <- arm_b_successes/nb
  WaldScore <- (pa_est-pb_est)/sqrt(pa_est*(1-pa_est)/(na-1) + pb_est*(1-pb_est)/(nb-1))
  reward <- (na*pa+nb*pb)/(na+nb)
  if(ap_output){
    da <- rbeta(10000,arm_a_successes[i-1]+1,arm_a_failures[i-1]+1)
    db <- rbeta(10000,arm_b_successes[i-1]+1,arm_b_failures[i-1]+1)
    pts <- mean(da>db)
    pc <- mean( abs(da-db)<c)
    pab <- pc/2 + (1-pc)*(pts)
    return(list(WaldScore=WaldScore,reward=reward,pab=pab,action_each_step=action_each_step,reward_each_step=reward_each_step,count=cbind(arm_a_successes,arm_a_failures,arm_b_successes,arm_b_failures)))
  } else{
    col_name <- c('WaldScore','reward','step','pop_p')
    output_df <- data.frame(array(dim=c(n,length(col_name))))
    colnames(output_df) <- col_name
    output_df['WaldScore'] <- WaldScore
    output_df['reward'] <- reward
    output_df['step'] <- 1:n
    output_df['pop_p'] <- (arm_a_successes+arm_b_successes)/(1:n)
    return(output_df)
    
    #wald, reward, population p
    #return(list(WaldScore=WaldScore,reward=((na*pa+nb*pb)/n),pop_p=(arm_a_s+arm_b_s)/n))
    #return(list(WaldScore=WaldScore,reward=reward,action_each_step=action_each_step,reward_each_step=reward_each_step,count=cbind(arm_a_successes,arm_a_failures,arm_b_successes,arm_b_failures)))
  }
  #return(c(WaldScore,arm_a_successes,na,arm_b_successes,nb,count_PC))
  #return(list(WaldScore=WaldScore,na=na,nb=nb,sa=arm_a_successes,sb=arm_b_successes))
  #return(list(WaldScore=WaldScore,reward=reward,pab=pab,count=cbind(arm_a_successes,arm_a_failures,arm_b_successes,arm_b_failures)))
  #return(c(pab,WaldScore[n],reward[n]))
}

#a fcuntion from output to FPR/Power/Reward/ratio

FPR_thre <- function(n,c,eps,burnin=0,bb=5000,p_list=seq(0.3,0.7,0.05)){
  #FPR_thre(n,c,eps,bb=50)
  wald <- array(dim=c(n,length(p_list),bb))
  for (i in 1:length(p_list)){
    p <- p_list[i]
    for (j in 1:bb){
      res <- TS2b(pa=p,pb=p,n=n,c=c,eps=eps,burnin=burnin)
      wald[,i,j] <- res$WaldScore
    }
  }
 
  return(list(thre=apply(abs(wald),c(1,2),quantile,probs=0.95,na.rm=T),fpr=apply(abs(wald)>1.96,c(1,2),mean,na.rm=T)))
}




#sim 2 has histroy tracking and lambda, can cahnge back
sim3 <- function(para,B,bb=1000){
  #simulation for 2-arm settings  
  #the output action/reward history is time_step(row)*simulation trials(column)
  pa <- para$pa
  pb <- para$pb
  n <- para$n
  c <- para$c
  eps <- para$eps
  
  reward_df <- array(dim=c(B,n))
  wald_df <- array(dim=c(B,n))
  pop_p_df <- array(dim=c(B,n))
  #arm_count_df <- array(dim=c(B,4))
  
  #reward_df <- c()
  #wald_df <- c()
  #pop_p_df <- c()
  
  #action and reward history is full info. I save those to Google drive. 
  #action_history <- array(dim=c(B,n)) 
  #reward_history <- array(dim=c(B,n))
  #objective_scores <- matrix(array(dim=c(n,length(lambdas))),ncol=length(lambdas))
  #colnames(objective_scores) <- lambdas
  if (is.null(para$burnin)){
    burnin <- 1
  }else{
    burnin <- para$burnin
  }
  
  # Initializes the progress bar
  probar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                             total = B,
                             complete = "=",   # Completion bar character
                             incomplete = "-", # Incomplete bar character
                             current = ">",    # Current bar character
                             clear = FALSE,    # If TRUE, clears the bar when finish
                             width = 100)      # Width of the progress bar
  for(i in 1:B) {
    probar$tick()
    res <- TS2b(pa,pb,n,c,eps,burnin=burnin)
    reward_df[i,] <- res$reward
    wald_df[i,] <- res$WaldScore
    pop_p_df[i,] <- res$pop_p
    #arm_count_df[i,,] <- res$count
    #action_history[i,] <- res$action_each_step
    #reward_history[i,] <- res$reward_each_step
  }
  reward <- apply(reward_df,2,mean,na.rm=T)
  #reward <- mean(reward_df)
  
  
  #get wald thre
  pmin <- 0.2
  pmax <- 0.8
  p_list <- seq(pmin,pmax,(pmax-pmin)/10)
  
  res <- FPR_thre(n,c,eps,p_list=p_list,bb=bb)
  wald_thre <- res$thre
  
  tests <- array(dim=c(B,n))
  nw <- dim(wald_thre)[2]
  seq2 <- seq(0,0.9,0.1)
  con_mat <- cbind(1-seq2,seq2)
  for (k in 1:n){
    
    wald_mat <- rbind(wald_thre[k,-nw],wald_thre[k,-1])
    
    walds <- c(as.vector(con_mat%*%wald_mat),wald_thre[k,nw])
    pop_ps <- c(as.vector(con_mat%*% (rbind(p_list[-nw],p_list[-1])) ),pmax)
    
    
    for (i in 1:B){
      tests[i,k] <- abs(wald_df[i,k])>walds[which(pop_p_df[i,k]<=pop_ps)[1]]
    }
  }
  
  
  
  
  #prop_x1 <- t(t(apply(arm_count_df[,,c(1,2)],c(1,2),sum))/c(1:n))
  #success_x1 <- arm_count_df[,,1]
  #success_x2 <- arm_count_df[,,3]
  
  #prop_x1 <- apply(prop_x1 , 2,mean)
  return(list(reward=reward,power= apply(tests,2,mean,na.rm=T) ))
  
}

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

#work flow































#paste back ..... 

###########################################################################
###############################  garbage  #################################
###########################################################################
n <- 300
c <- 0.05
eps <- 0
pa <- 0.5
pb <- 0.5
fpr <- c()
for ( i in 1:1000){
  res <- TS2b(pa,pb,n,c,eps,burnin=burnin)
  fpr[i] <- abs(res[1])>FPR_thre(res[2],res[3],n,c,eps,bb=100)
}



sim <- function(para,B){
  #simulation for 2-arm settings  
  #the output action/reward history is time_step(row)*simulation trials(column)
  pa <- para$pa
  pb <- para$pb
  n <- para$n
  c <- para$c
  eps <- para$eps
  
  probar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                             total = B,
                             complete = "=",   # Completion bar character
                             incomplete = "-", # Incomplete bar character
                             current = ">",    # Current bar character
                             clear = FALSE,    # If TRUE, clears the bar when finish
                             width = 100)      # Width of the progress bar
  for(i in 1:B) {
    probar$tick()
    res <- TS2b(pa,pb,n,c,eps,burnin=burnin)
    reward_df[i,] <- res$reward
    wald_df[i,] <- res$WaldScore
    arm_count_df[i,,] <- res$count
    action_history[i,] <- res$action_each_step
    reward_history[i,] <- res$reward_each_step
  }
  reward <- apply(reward_df,2,mean,na.rm=T)
  if (abs(pa-pb)>0){
    if(pa>pb){
      wald_test <- apply(wald_df>1.96,2,mean,na.rm=T)
    } else{
      wald_test <- apply(wald_df< -1.96,2,mean,na.rm=T)
    }
  }else{
    wald_test <- apply(abs(wald_df)>1.96,2,mean,na.rm=T)
  }
  prop_x1 <- t(t(apply(arm_count_df[,,c(1,2)],c(1,2),sum))/c(1:n))
  success_x1 <- arm_count_df[,,1]
  success_x2 <- arm_count_df[,,3]
  for (i in 1:length(lambdas)){
    ress <- entropy_objective_optimal_2arm(lambda=lambdas[i],delta=pa-pb,prop_x1)
    objective_scores[,i] <- apply(ress$loss,2,mean,na.rm=T)
  }
  
}

sim2 <- function(para,B,lambdas=NA){
#simulation for 2-arm settings  
#the output action/reward history is time_step(row)*simulation trials(column)
  pa <- para$pa
  pb <- para$pb
  n <- para$n
  c <- para$c
  eps <- para$eps
  
  reward_df <- array(dim=c(B,n))
  wald_df <- array(dim=c(B,n))
  arm_count_df <- array(dim=c(B,n,4))
  
  #action and reward history is full info. I save those to Google drive. 
  action_history <- array(dim=c(B,n)) 
  reward_history <- array(dim=c(B,n))
  #objective_scores <- matrix(array(dim=c(n,length(lambdas))),ncol=length(lambdas))
  #colnames(objective_scores) <- lambdas
  if (is.null(para$burnin)){
    burnin <- 1
  }else{
    burnin <- para$burnin
  }

  # Initializes the progress bar
  probar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = B,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  for(i in 1:B) {
    probar$tick()
    res <- TS2b(pa,pb,n,c,eps,burnin=burnin)
    reward_df[i,] <- res$reward
    wald_df[i,] <- res$WaldScore
    arm_count_df[i,,] <- res$count
    action_history[i,] <- res$action_each_step
    reward_history[i,] <- res$reward_each_step
  }
  reward <- apply(reward_df,2,mean,na.rm=T)
  if (abs(pa-pb)>0){
    if(pa>pb){
      wald_test <- apply(wald_df>1.96,2,mean,na.rm=T)
    } else{
      wald_test <- apply(wald_df< -1.96,2,mean,na.rm=T)
    }
  }else{
    wald_test <- apply(abs(wald_df)>1.96,2,mean,na.rm=T)
  }
  prop_x1 <- t(t(apply(arm_count_df[,,c(1,2)],c(1,2),sum))/c(1:n))
  success_x1 <- arm_count_df[,,1]
  success_x2 <- arm_count_df[,,3]
  for (i in 1:length(lambdas)){
    ress <- entropy_objective_optimal_2arm(lambda=lambdas[i],delta=pa-pb,prop_x1)
    objective_scores[,i] <- apply(ress$loss,2,mean,na.rm=T)
  }
  #prop_x1 <- apply(prop_x1 , 2,mean)
  return(list(reward=reward,power=wald_test,prop_x1 =prop_x1,
              s1=success_x1,s2=success_x2,
              #objective_scores=objective_scores,
              action_history=as.data.frame(t(action_history)),
              reward_history=as.data.frame(t(reward_history))))
  
}











