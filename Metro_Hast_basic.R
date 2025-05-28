# Function to calculate the value of the log likelihood
#Function to simulate data 
make.data <- function(lambda,rho,t){
  a2=sample(seq(1,t,1), 400, replace=TRUE)
  young_prob <- 1-exp(-rho*lambda*a2)
  a1=sample(seq(t+1,70,1), 400, replace=TRUE)
  old_prob <- 1-exp(lambda*(t-a1)-t*rho*lambda)
  a <- c(a2,a1)
  prob <- c(young_prob,old_prob)
  
  x=vector()
  for(i in 1:length(prob)){
    x[i]=rbinom(1,1,prob[i])
  }
  return(list(a,x))
}


rep=1
sim <- data.frame(lambda=sample(seq(0.001,0.005,0.0001),rep,replace=TRUE), 
                  rho=sample(seq(2,6,0.001),rep,replace=TRUE), 
                  t=sample(seq(10,30,0.1),rep,replace=TRUE),
                  lambda_hat = rep(NA,rep), rho_hat=rep(NA,rep), t_hat=rep(NA,rep),
                  lambda_hat2 = rep(NA,rep), rho_hat2=rep(NA,rep), t_hat2=rep(NA,rep))

d <- with(sim, make.data(lambda[1], rho[1], t[1]))
  a <- d[[1]]
  x <- d[[2]]
  

lambda_start=0.0025
rho_start=4
t_start=20
iter=5000
a
x
# Function to run the Metropolis-Hastings Algorithm

Metro_Hast <- function (lambda_start,rho_start,t_start,iter,a,x) {
  lambda=rep(NA, iter) 
  lambda[1]=lambda_start
  rho=rep(NA, iter)
  rho[1]=rho_start
  t=rep(NA, iter)
  t[1]=t_start
  loglik=rep(NA, iter)
  loglik[1]=if(lambda[1]<0|lambda[1]>0.1|rho[1]<1|rho[1]>15|t[1]<5|t[1]>65){
    -Inf
  } else{
    liknum=vector()
    for(j in 1:length(a)){
      if(t[1]<a[j]){
        liknum[j]=log(exp((lambda[1]*(t[1]-a[j])-t[1]*rho[1]*lambda[1])*abs(x[j]-1))*(1-exp(lambda[1]*(t[1]-a[j])-t[1]*rho[1]*lambda[1]))^x[j])
      }else{
        liknum[j]=log(exp(-rho[1]*lambda[1]*a[j]*abs(x[j]-1))*(1-exp(-rho[1]*lambda[1]*a[j]))^x[j])
      }
    }
    sum(liknum)
  }
  accept=rep(NA, iter)
  accept[1]=0
  
  prop.sd <- c(0.0001,0.5,10)
  
  for (i in 2:iter) {
    mu <- c(lambda[i-1], rho[i-1], t[i-1])
    theta.new <- rnorm(c(1,1,1),mu, prop.sd)
      
    loglik.new <- if(theta.new[1]<0|theta.new[1]>0.1|theta.new[2]<1|theta.new[2]>15|theta.new[3]<5|theta.new[3]>65){
        -Inf
      } else{
        liknum=vector()
        for(j in 1:length(a)){
          if(theta.new[3]<a[j]){
            liknum[j]=log(exp((theta.new[1]*(theta.new[3]-a[j])-theta.new[3]*theta.new[2]*theta.new[1])*abs(x[j]-1))*(1-exp(theta.new[1]*(theta.new[3]-a[j])-theta.new[3]*theta.new[2]*theta.new[1]))^x[j])
          }else{
            liknum[j]=log(exp(-theta.new[2]*theta.new[1]*a[j]*abs(x[j]-1))*(1-exp(-theta.new[2]*theta.new[1]*a[j]))^x[j])
          }
        }
        sum(liknum)
      }
      
    log.ratio <- loglik.new - loglik[i-1]
    log.accept.prob <- min(log.ratio, 0)
    
    if (log (runif (1)) < log.accept.prob) {
      lambda[i]=theta.new[1]
      rho[i]=theta.new[2]
      t[i]=theta.new[3]
      loglik[i]=loglik.new
      accept[i]=1
    } else {
      lambda[i]=lambda[i-1]
      rho[i]=rho[i-1]
      t[i]=t[i-1]
      loglik[i]=loglik[i-1]
      accept[i]=0
    } 
  }
  return(data.frame(lambda,rho,t, loglik, accept))
}

start_time <- Sys.time()
try <- Metro_Hast(0.0025, 4, 20,5000, a,x)
end_time <- Sys.time()
end_time - start_time 

library(Rcpp)

sourceCpp("Metro_Hast_basic.cpp")
sourceCpp("calcLogLik.cpp")

start_time <- Sys.time()
try <- Metro_Hast(0.0025, 4, 20,5000, a,x)
end_time <- Sys.time()
end_time - start_time 


