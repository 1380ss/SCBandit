#Goal: givne output, 
#load that function TS2b
#1. function calculating optimization of objective with differnt lambda
# look at your overleaf..
#2. calculate empirical performance (at each step, each 'c' and epsilon.)

source("Algorithm for two arm.R")

#function calculating optimal solution
entropy_objective_optimal_2arm <- function(lambda,delta,x=NaN,min_sacle=0){
  #delta = mu1-mu2, x is proportion to arm 1
  #objective = delta x + lambda [xlogx + (1-x) log(1-x)]
  xx <- exp(delta/lambda)
  if(is.na(x)){
    x <- xx/(1+xx)
  }
  scale_factor <- (lambda*log(1+xx))
  if(is.na(min_sacle)){
    min_sacle <- delta*0.5+lambda*log(2)
  }
  
  loss <- (x*delta-lambda*( x*log(x)+(1-x)*log(1-x))-min_sacle)/(scale_factor-min_sacle)
  return(list(best_x=xx/(1+xx),loss=loss))
}
entropy_objective_optimal_2arm(lambda=0.2,0.2,x=0.9)
