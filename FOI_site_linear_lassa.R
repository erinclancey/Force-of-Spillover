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
Sites=unique(sero_data$district)
# Estimate FOI (lambda) by site
df <- create_empty_table(length(Sites), 10)
colnames(df) <-  c("lambda","lambdaSE","beta","betaSE","t","tSE","loglik", "Site","Size","Seroprev")


for(i in 1:length(Sites)){
  site=subset(sero_data, district==Sites[i])
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
      beta=theta[2]
      tau=theta[3]
      if (lambda<0) {
        return(Inf)
      } else {
        A<-max(a)
        array=cbind(a,x)
        old=array[tau-(A-a)>0,]
        young=array[tau-(A-a)<=0,]
        lik_old=sum( log (
          exp((-old[,1]*lambda/1000-1/2*beta/10000*(A-tau)^2)*(1-old[,2]))*
            (1-exp(-old[,1]*lambda/1000-1/2*beta/10000*(A-tau)^2))^old[,2] ) )
        lik_young=sum( log( 
          exp(young[,1]*(beta/10000*(young[,1]/2-A-tau)-lambda/1000)*(1-young[,2]))*
            (1-exp(young[,1]*(beta/10000*(young[,1]/2-A-tau)-lambda/1000)))^young[,2] ) )
        
        loglik=lik_old+lik_young
        return(-loglik)
      }
    }
    #Maximize the Log Likelihood for Lambda (L)
    m <- optim(c(0.1,0.001,70), nloglik, method = "SANN", hessian = TRUE)
    df[i,1]=m$par[1]/1000
    df[i,2]=sqrt(diag(solve(m$hessian)))[1]/1000
    df[i,3]=m$par[2]/10000
    df[i,4]=sqrt(diag(solve(m$hessian)))[2]/10000
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

