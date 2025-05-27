library(tidyverse)
library(MASS)
library(Boom)
library(coda)
library(bayestestR)
library(lhs)
library(MCMCglmm)
library(pracma)
library(latex2exp)
library(ggpubr)
library(optimx)
#Function to simulate data 
make.data <- function(lambda,beta,tau){
  a <- sample(seq(1,80,0.1), 1000, replace=TRUE)
  A=max(a)
  prob <- vector()
  for(i in 1:length(a)){
    if(tau-(A-a[i])>0){
      prob[i] <- 1-exp(-a[i]*lambda-1/2*beta*(A-tau)^2)
    }else{
      prob[i] <- 1-exp(a[i]*(beta*(a[i]/2-A-tau)-lambda))
    }
  }
  x=vector()
  for(i in 1:length(prob)){
    x[i]=rbinom(1,1,prob[i])
  }
  return(list(a,x,prob))
}

nloglik <- function(theta){
  lambda=theta[1]
  beta=theta[2]
  tau=theta[3]
  if (lambda<0) {
    return(Inf)
  } else {
    lik_old=sum( log (
        exp((-old[,1]*lambda-1/2*beta*(A-tau)^2)*(1-old[,2]))*
          (1-exp(-old[,1]*lambda-1/2*beta*(A-tau)^2))^old[,2] ) )
    
    lik_young=sum( log( 
        exp(young[,1]*(beta*(young[,1]/2-A-tau)-lambda)*(1-young[,2]))*
          (1-exp(young[,1]*(beta*(young[,1]/2-A-tau)-lambda)))^young[,2] ) )
  
    loglik=lik_old+lik_young
    return(-loglik)
  }
}

rep=1
sim <- data.frame(lambda=sample(seq(0.001,0.005,0.0001),rep,replace=TRUE), 
                  beta=sample(seq(0.0001, 0.0005,0.00001),rep,replace=TRUE), 
                  tau=sample(seq(40,50,0.1),rep,replace=TRUE),
                  lambda_hat = rep(NA,rep), beta_hat=rep(NA,rep), tau_hat=rep(NA,rep),
                  seroprev=rep(NA,rep),converge=rep(NA,rep), maxit=rep(NA,rep))

for(j in 1:nrow(sim)){
  
  d <- with(sim, make.data(lambda[j], beta[j], tau[j]))
  
  a <- d[[1]]
  x <- d[[2]]
  A<-max(a)
  array=cbind(a,x)
  old=array[tau-(A-a)>0,]
  young=array[tau-(A-a)<=0,]
  #young <- matrix(young,ncol=2)
  
  sim$seroprev[j] <- sum(x)/length(x)
  start <- c(0.001, 0.0001, 40)
  
  m <- optimx(start, nloglik, method = "nlminb",
              lower = c(0.0001, 0.00001, 30), upper = c(0.1, 0.01, 70),
              control =list(iter.max = 1000, eval.max = 1000, trace=0,
                            xf.tol=2.2e-20,rel.tol=1e-5,step.min=1, step.max=1))
  
  sim$lambda_hat[j] <- m$p1
  sim$beta_hat[j] <- m$p2
  sim$tau_hat[j] <- m$p3
  sim$converge[j] <- m$convcode
  sim$maxit[j] <- m$niter
  
  print(paste(round(j/nrow(sim) * 100, 1), "%", sep = ""))
}

p1 <- ggplot(sim, aes(x=lambda/1000, y=lambda_hat/1000)) + 
  geom_point(color="#A73030FF", shape=20, size=3) + 
  xlab(TeX("$\\lambda$")) + ylab(TeX("$\\hat{\\lambda}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=1, color="black")+
  xlim(c(1e-04,0.008))+ylim(c(1e-04,0.008))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p2 <- ggplot(sim, aes(x=rho, y=rho_hat)) + 
  geom_point( color="#0073C2FF",shape=20, size=3) + 
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=1, color="black")+
  xlim(c(1,10))+ylim(c(1,10))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p3 <- ggplot(sim, aes(x=t, y=t_hat)) +
  geom_point( color="#868686FF",shape=20, size=3) + 
  xlab(TeX("$t$")) + ylab(TeX("$\\hat{t}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=1, color="black")+
  xlim(c(5,25))+ylim(c(5,25))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
ggarrange(p1,p2,p3, 
          ncol = 3, nrow = 1,widths=c(2.2,2,2))