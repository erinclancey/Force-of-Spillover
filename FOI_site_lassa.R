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
df <- create_empty_table(length(Sites), 6)
colnames(df) <-  c("lambda","se","loglik", "Site","Size","Serprev")


for(i in 1:length(Sites)){
  site=subset(sero_data, chiefdom==Sites[i])
  age=site$age_en
  result=site$iggbinary
  if(sum(result)==0){
    df[i,1]=0
    df[i,2]=NA
    df[i,3]=NA
    df[i,4]=Sites[i]
  }else{
    nloglik <- function(lambda) {    
      if (lambda<0) {
        return(Inf)
      } else {
        loglik <- sum(log(exp(-lambda/100*age*(1-result))*(1-exp(-lambda/100*age))^result))
        return(-loglik)
      }
    }
    #Maximize the Log Likelihood for Lambda (L)
    m <- optim(c(lambda=0.1), nloglik, method = "SANN", hessian = TRUE)
    df[i,1]=m$par/100
    df[i,2]=sqrt(diag(solve(m$hessian)))/100
    df[i,3]=m$value
    df[i,4]=Sites[i]
  }
  df[i,5]=nrow(site)
  df[i,6]=sum(result)/nrow(site)
}
#view results by site
df
