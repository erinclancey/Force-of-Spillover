setwd("~/FOS Lassa")

library(tidyverse)
sero_data <- read.csv("Lassa_sero.csv", header=TRUE)
sero_data <- na.omit(sero_data)
create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}




#Summarize Serology by site and age
SiteByAge=sero_data %>%
  group_by(chiefdom) %>%
  summarize(mean_RVFserpos = mean(iggbinary, na.rm = TRUE),mean_Age = mean(age_en, na.rm = TRUE))
#Check for relationship between average age and IgG
ggplot(SiteByAge, aes(x=mean_Age, y=mean_RVFserpos)) + geom_point()

# Estimate FOI (lambda) by site
Sites=unique(sero_data$chiefdom)
# Estimate FOI (lambda) by site
df <- create_empty_table(length(Sites), 10)
colnames(df) <-  c("lambda","lambdaSE","rho","rhoSE","t","tSE","loglik", "Site","Size","Seroprev")


for(i in 1:length(Sites)){
  site=subset(sero_data, chiefdom==Sites[i])
  a=site$age_en
  x=site$iggbinary
  array=cbind(a,x)
  
  if(sum(x)==0){
    df[i,1]=0
    df[i,2]=NA
    df[i,3]=NA
    df[i,4]=NA
    df[i,5]=NA
    df[i,6]=NA
    df[i,7]=NA
    df[i,8]=Sites[i]
  }else{
    nloglik <- function(theta){
      lambda=theta[1]
      rho=theta[2]
      t=theta[3]
      if (lambda<0|rho<0|t<1) {
        return(Inf)
      } else {
      old=array[t<a,]
      young=array[t>=a,]
      young <- matrix(young,ncol=2)
          lik_old=sum(log(exp((lambda/100*(t-old[,1])-t*rho*lambda/100)*(1-old[,2]))*(1-exp(lambda/100*(t-old[,1])-t*rho*lambda/100))^old[,2]))
          lik_young=sum(log(exp(-rho*lambda/100*young[,1]*(1-young[,2]))*(1-exp(-rho*lambda/100*young[,1]))^young[,2]))
          loglik=lik_old+lik_young
      return(-loglik)
      }
    }
    #Maximize the Log Likelihood for Lambda (L)
    m <- optim(c(0.1,1,10), nloglik, method = "SANN", hessian = TRUE)
    df[i,1]=m$par[1]/100
    df[i,2]=sqrt(diag(solve(m$hessian)))[1]/100
    df[i,3]=m$par[2]
    df[i,4]=sqrt(diag(solve(m$hessian)))[2]
    df[i,5]=m$par[3]
    df[i,6]=sqrt(diag(solve(m$hessian)))[3]
    df[i,7]=m$value
    df[i,8]=Sites[i]
  }
  df[i,9]=nrow(site)
  df[i,10]=sum(x)/nrow(site)
}

#view results by site
df

