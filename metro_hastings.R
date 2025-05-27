
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

# Function to calculate the value of the log likelihood
calcLogLik <- function(lambda,rho,t,a,x){
  if (lambda<0|lambda>0.1|rho<1|rho>15|t<5|t>65){
    loglik <- -Inf
  } else{
    lik=vector()
    for(i in 1:length(a)){
      if(t<a[i]){
        lik[i]=log(exp((lambda*(t-a[i])-t*rho*lambda)*abs(x[i]-1))*(1-exp(lambda*(t-a[i])-t*rho*lambda))^x[i])
      }else{
        lik[i]=log(exp(-rho*lambda*a[i]*abs(x[i]-1))*(1-exp(-rho*lambda*a[i]))^x[i])
      }
    }
    loglik=sum(lik)
  }
  return(loglik)
}

# Function to run the Metropolis-Hastings Algorithm
get_proposal <- function (lambda_start,rho_start,t_start,iter,a,x) {
  lambda=rep(NA, iter) 
  lambda[1]=lambda_start
  rho=rep(NA, iter)
  rho[1]=rho_start
  t=rep(NA, iter)
  t[1]=t_start
  lik=rep(NA, iter)
  lik[1]=calcLogLik(lambda_start,rho_start,t_start, a,x)
  accept=rep(NA, iter)
  accept[1]=0
  
  prop.sd <- c(0.0001,0.5,10)
  
  for (i in 2:iter) {
    mu <- c(lambda[i-1], rho[i-1], t[i-1])
    theta.new <- rnorm(c(1,1,1),mu, prop.sd)
    loglik.new <- calcLogLik(theta.new[1],theta.new[2],theta.new[3], a,x)
    loglik.old <- calcLogLik(lambda[i-1], rho[i-1], t[i-1], a,x)
    log.ratio <- loglik.new - loglik.old
    log.accept.prob <- min(log.ratio, 0)
    
    if (log (runif (1)) < log.accept.prob) {
      lambda[i]=theta.new[1]
      rho[i]=theta.new[2]
      t[i]=theta.new[3]
      lik[i]=loglik.new
      accept[i]=1
    } else {
      lambda[i]=lambda[i-1]
      rho[i]=rho[i-1]
      t[i]=t[i-1]
      lik[i]=loglik.old
      accept[i]=0
    } 
  }
  p_cur=sum(accept)/iter
  prop.sd_new <- cov(cbind(lambda,rho,t))*qnorm(0.238/2, mean=0, sd=1)/qnorm(p_cur/2)
  return(prop.sd_new)
}

Metro_Hast2 <- function (lambda_start,rho_start,t_start,iter1,iter2,a,x) {
  lambda=rep(NA, iter) 
  lambda[1]=lambda_start
  rho=rep(NA, iter)
  rho[1]=rho_start
  t=rep(NA, iter)
  t[1]=t_start
  lik=rep(NA, iter)
  lik[1]=calcLogLik(lambda_start,rho_start,t_start, a,x)
  accept=rep(NA, iter)
  accept[1]=0
  
  prop.sd <- prop.sd_new
  
  for (i in 2:iter1) {
    mu <- c(lambda[i-1], rho[i-1], t[i-1])
    theta.new <- mvrnorm(1,mu, prop.sd)
    loglik.new <- calcLogLik(theta.new[1],theta.new[2],theta.new[3], a,x)
    loglik.old <- calcLogLik(lambda[i-1], rho[i-1], t[i-1], a,x)
    log.ratio <- loglik.new - loglik.old
    log.accept.prob <- min(log.ratio, 0)
    
    if (log (runif (1)) < log.accept.prob) {
      lambda[i]=theta.new[1]
      rho[i]=theta.new[2]
      t[i]=theta.new[3]
      lik[i]=loglik.new
      accept[i]=1
    } else {
      lambda[i]=lambda[i-1]
      rho[i]=rho[i-1]
      t[i]=t[i-1]
      lik[i]=loglik.old
      accept[i]=0
    } 
  }

  
  return(list(lambda,rho,t, lik, accept))
}

############Start simulations
                  
rep=1
sim <- data.frame(lambda=sample(seq(0.001,0.005,0.0001),rep,replace=TRUE), 
                  rho=sample(seq(2,6,0.001),rep,replace=TRUE), 
                  t=sample(seq(10,30,0.1),rep,replace=TRUE),
                  lambda_hat = rep(NA,rep), rho_hat=rep(NA,rep), t_hat=rep(NA,rep),
                  lambda_hat2 = rep(NA,rep), rho_hat2=rep(NA,rep), t_hat2=rep(NA,rep))
start_time <- Sys.time()
for(j in 1:nrow(sim)){
  
  d <- with(sim, make.data(lambda[j], rho[j], t[j]))
  
  a <- d[[1]]
  x <- d[[2]]
  
  get_proposal(0.0025, 4, 20,5000, a,x)
  
  results1 <- Metro_Hast1(0.0025, 4, 20, iter=5000, a,x)
  results.post1 <- do.call(cbind, results1)
  post1 <- data.frame(results.post1)
  colnames(post1) <- c('lambda','rho','t','loglik','accept')
  post_theta1 <- post1[,1:3]
  p_cur=sum(post1$accept)/1999
  prop.sd_new <- cov(post_theta1)*qnorm(0.238/2, mean=0, sd=1)/qnorm(p_cur/2)
  
  #post.filtered = post1[seq(1, nrow(post1), 30), ]
  
  
  results2 <- Metro_Hast2(mean(post1$lambda), mean(post1$rho), mean(post1$t), iter=10000, a,x)
  results.post2 <- do.call(cbind, results2)
  post2 <- data.frame(results.post2)
  colnames(post2) <- c('lambda','rho','t','loglik','accept')
  post.filtered = post2[seq(1, nrow(post2), 50), ]
  
  sim$lambda_hat[j] <- Mode(post.filtered$lambda)
  sim$rho_hat[j] <- Mode(post.filtered$rho)
  sim$t_hat[j] <- Mode(post.filtered$t)
  
  sim$lambda_hat2[j] <- mean(post.filtered$lambda)
  sim$rho_hat2[j] <- mean(post.filtered$rho)
  sim$t_hat2[j] <- mean(post.filtered$t)
  
  print(paste(round(j/nrow(sim) * 100, 1), "%", sep = ""))
}
end_time <- Sys.time()
end_time - start_time 

p1 <- ggplot(sim, aes(x=lambda, y=lambda_hat)) +
  geom_point(color="#A73030FF", shape=20, size=3) +
  xlab(TeX("$\\lambda$")) + ylab(TeX("$\\hat{\\lambda}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#0073C2FF")+
  theme_bw()
p2 <- ggplot(sim, aes(x=rho, y=rho_hat)) +
  geom_point( color="#EFC000FF",shape=20, size=3) +
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#0073C2FF")+
  theme_bw()
p3 <- ggplot(sim, aes(x=t, y=t_hat)) +
  geom_point( color="#868686FF",shape=20, size=3) +
  xlab(TeX("$t$")) + ylab(TeX("$\\hat{t}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="#0073C2FF")+
  theme_bw()
ggarrange(p1,p2,p3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


p1 <- ggplot(sim, aes(x=lambda, y=lambda_hat2)) + 
  geom_point(color="#A73030FF", shape=20, size=3) + 
  xlab(TeX("$\\lambda$")) + ylab(TeX("$\\hat{\\lambda}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=1, color="#0073C2FF")+
  ylim(c(0,0.01))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p2 <- ggplot(sim, aes(x=rho, y=rho_hat2)) + 
  geom_point( color="#EFC000FF",shape=20, size=3) + 
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=1, color="#0073C2FF")+
  ylim(c(0,15))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p3 <- ggplot(sim, aes(x=t, y=t_hat2)) +
  geom_point( color="#868686FF",shape=20, size=3) + 
  xlab(TeX("$t$")) + ylab(TeX("$\\hat{t}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=1, color="#0073C2FF")+
  ylim(c(5,40))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
ggarrange(p1,p2,p3, 
          ncol = 3, nrow = 1,widths=c(2.2,2,2))

"#A73030FF"

plot(post1$lambda, type="l")
plot(post1$rho, type="l")
plot(post1$t, type="l")

hist(post1$lambda)
hist(post1$rho)
hist(post1$t)

acf(post1$lambda)
acf(post1$rho)
acf(post1$t)

plot(post.filtered$lambda, type="l")
plot(post.filtered$rho, type="l")
plot(post.filtered$t, type="l")

hist(post.filtered$lambda)
hist(post.filtered$rho)
hist(post.filtered$t)

acf(post.filtered$lambda)
acf(post.filtered$rho)
acf(post.filtered$t)

ggplot(post.filtered, aes(x=lambda, y=rho))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75, adjust=2)+ 
  geom_density_2d(colour="black", adjust=2)+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste(lambda))) +  
  ylab(expression(paste(rho))) 

ggplot(post.filtered, aes(x=lambda, y=t))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75, adjust=2)+ 
  geom_density_2d(colour="black", adjust=2)+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste(lambda))) +  
  ylab(expression(paste(t))) 

ggplot(post.filtered, aes(x=rho, y=t))+
  geom_density_2d_filled(show.legend = TRUE, alpha=0.75, adjust=2)+ 
  geom_density_2d(colour="black", adjust=2)+
  geom_point(alpha = 0,size=0.5, show.legend=FALSE) +theme_minimal()+
  xlab(expression(paste(rho))) +  
  ylab(expression(paste(t)))





