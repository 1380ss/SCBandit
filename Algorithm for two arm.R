#library lists

#no I shouldn't have library loaded here...

#TS algorithm for 2-arm, binary reward case. Can use PostDiff and epsilon-UR
TS2b <- function(pa,pb,n,c,eps,burnin=1,ap_output=F){
  # n is the sample size (785 or 197 in our case)
  arm_a_successes <- c(0)
  arm_a_failures <- c(0)
  arm_b_successes <- c(0)
  arm_b_failures <- c(0)
  
  arm_a_s <- c(0)
  arm_a_f <- c(0)
  arm_b_s <- c(0)
  arm_b_f <- c(0)
  
  for (i in 1:burnin){
    draw_a <- runif(1)
    draw_b <- runif(1)
    if (draw_a>draw_b){
      if(runif(1)<pa){
        arm_a_s <- arm_a_s+1
      } else {
        arm_a_f <- arm_a_f+1
      }
    } else{
      if(runif(1)<pb){
        arm_b_s <- arm_b_s+1
      } else {
        arm_b_f <- arm_b_f+1
      }
    }
    arm_a_successes[i] <- arm_a_s
    arm_a_failures[i] <- arm_a_f
    arm_b_successes[i] <- arm_b_s
    arm_b_failures[i] <- arm_b_f
  }
  
  for (i in (burnin+1):n){
    if (runif(1)<eps){
      #epsilon-UR step
      draw_a <- runif(1)
      draw_b <- runif(1)
    } else {
      #TS-step
      draw_a <- rbeta(1,arm_a_successes[i-1]+1,arm_a_failures[i-1]+1)
      draw_b <- rbeta(1,arm_b_successes[i-1]+1,arm_b_failures[i-1]+1)
      if (abs(draw_a-draw_b)<c){
        #PostDiff step
        draw_a <- runif(1)
        draw_b <- runif(1)
      } else {
        draw_a <- rbeta(1,arm_a_successes[i-1]+1,arm_a_failures[i-1]+1)
        draw_b <- rbeta(1,arm_b_successes[i-1]+1,arm_b_failures[i-1]+1)
      }
    }
    if (draw_a>draw_b){
      if(runif(1)<pa){
        arm_a_s <- arm_a_s+1
      } else {
        arm_a_f <- arm_a_f+1
      }
    } else{
      if(runif(1)<pb){
        arm_b_s <- arm_b_s+1
      } else {
        arm_b_f <- arm_b_f+1
      }
    }
    arm_a_successes[i] <- arm_a_s
    arm_a_failures[i] <- arm_a_f
    arm_b_successes[i] <- arm_b_s
    arm_b_failures[i] <- arm_b_f
  }
  #na: number of times arm A was pulled. nb: same
  
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
    return(list(WaldScore=WaldScore,reward=reward,pab=pab,count=cbind(arm_a_successes,arm_a_failures,arm_b_successes,arm_b_failures)))
  } else{
    return(list(WaldScore=WaldScore,reward=reward,count=cbind(arm_a_successes,arm_a_failures,arm_b_successes,arm_b_failures)))
  }
  #return(c(WaldScore,arm_a_successes,na,arm_b_successes,nb,count_PC))
  #return(list(WaldScore=WaldScore,na=na,nb=nb,sa=arm_a_successes,sb=arm_b_successes))
  #return(list(WaldScore=WaldScore,reward=reward,pab=pab,count=cbind(arm_a_successes,arm_a_failures,arm_b_successes,arm_b_failures)))
  #return(c(pab,WaldScore[n],reward[n]))
}

#a fcuntion from output to FPR/Power/Reward/ratio
para <- list(pa=0.6,pb=0.4,n=785,c=0.1,eps=0)
sim_2b <- function(para,B){
  pa <- para$pa
  pb <- para$pb
  n <- para$n
  c <- para$c
  eps <- para$eps
  if (is.null(para$burnin)){
    burnin <- 1
  }else{
    burnin <- para$burnin
  }
  TS2b(pa,pb,n,c,eps,burnin=burnin)
}
