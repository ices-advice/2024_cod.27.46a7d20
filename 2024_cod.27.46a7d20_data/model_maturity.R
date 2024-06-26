library(lme4)      
library(DHARMa)
library(reshape2)  
library(ggplot2)   
library(xtable)   
library(patchwork) 
library(glmmTMB)

load("data/mat_data.Rdata")

working_year <- 2024
years <- 1983:working_year
ages <- sort(unique(mat_data$Age))

mo <- list()

#########################
# Maturity ogive - South
#########################

dat = subset(mat_data,SubArea=="S")

dat$wt = dat$wi*dat$wtC
dat$y = dat$NoMature
dat$n = dat$NoAtAlk
dat$fyear = factor(dat$Year)
dat$Year = as.numeric(as.character(dat$fyear))  
dat$fage=factor(dat$Age)
dat$Age = dat$Age


## M4 has a fixed age effects and iid random effects for year and age-year devs. Not really a model but a data summary

m4 <- glmmTMB(cbind(y,n-y) ~ fage - 1+ (1|fyear) + (1|(fage:fyear)),
              data  = dat, weight=wt,family = "binomial")

pred.data2=expand.grid(Age=ages,Year=years)  
pred.data2$fage = factor(pred.data2$Age)        
pred.data2$fyear = factor(pred.data2$Year)       
pred.data2$wt = 1

tmp = predict(m4,type="response",newdata=pred.data2,se.fit=TRUE,allow.new.levels=TRUE)
pred.data2$pred=tmp$fit
pred.data2$pred.se = tmp$se.fit  
tmp = predict(m4,type="link",newdata=pred.data2,se.fit=TRUE,allow.new.levels=TRUE)  
pred.data2$pred.logit=tmp$fit
pred.data2$pred.logit.se = tmp$se.fit

pred.data2$logit.L = pred.data2$pred.logit - 1.96*pred.data2$pred.logit.se
pred.data2$logit.U = pred.data2$pred.logit + 1.96*pred.data2$pred.logit.se  
pred.data2$L = exp(pred.data2$logit.L)/(1+exp(pred.data2$logit.L))      
pred.data2$U = exp(pred.data2$logit.U)/(1+exp(pred.data2$logit.U))
pred.data2$U[is.na(pred.data2$U)] =1 
pred.data2$logit.U[pred.data2$logit.U>40] =40
pred.data2$logit.L[pred.data2$logit.L< (-40)] =-40

ggplot(pred.data2,aes(x=Year,y=pred))+  
  facet_wrap(~Age, ncol = 2,scales="free_y") + 
  geom_errorbar(aes(ymin=L,ymax=U),width=0.2,color='green')+    
  geom_line()+  
  xlab("Year")+ylab("Proportion Mature, with CI") +
  theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
        legend.position="top",legend.key.width=unit(2, 'cm'),
        legend.text = element_text(size=8)) +
  theme_bw()

# ggplot(pred.data2,aes(x=Year,y=pred.logit))+
#   facet_wrap(~Age, ncol = 2,scales="free_y") +
#   geom_errorbar(aes(ymin=logit.L,ymax=logit.U),width=0.2,color='grey')+
#   geom_line()+
#   xlab("Year")+ylab("Logit Mature, with CI") +
#   theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
#         legend.position="top",legend.key.width=unit(2, 'cm'),
#         legend.text = element_text(size=8)) +
#   theme_bw()

mo[[1]] <- acast(pred.data2, fyear ~ fage, value.var = "pred")


############################
# Maturity ogive - Northwest
############################

rm(list=ls()[!ls() %in% c("mat_data","ages","years","mo")])

dat = subset(mat_data,SubArea=="NW")

dat$wt = dat$wi*dat$wtC
dat$y = dat$NoMature
dat$n = dat$NoAtAlk
dat$fyear = factor(dat$Year)
dat$Year = as.numeric(as.character(dat$fyear))  
dat$fage=factor(dat$Age)
dat$Age = dat$Age


## M4 has a fixed age effects and iid random effects for year and age-year devs. Not really a model but a data summary

m4 <- glmmTMB(cbind(y,n-y) ~ fage - 1+ (1|fyear) + (1|(fage:fyear)),
              data  = dat, weight=wt,family = "binomial")

pred.data2=expand.grid(Age=ages,Year=years)  
pred.data2$fage = factor(pred.data2$Age)        
pred.data2$fyear = factor(pred.data2$Year)       
pred.data2$wt = 1

tmp = predict(m4,type="response",newdata=pred.data2,se.fit=TRUE,allow.new.levels=TRUE)
pred.data2$pred=tmp$fit
pred.data2$pred.se = tmp$se.fit  
tmp = predict(m4,type="link",newdata=pred.data2,se.fit=TRUE,allow.new.levels=TRUE)  
pred.data2$pred.logit=tmp$fit
pred.data2$pred.logit.se = tmp$se.fit

pred.data2$logit.L = pred.data2$pred.logit - 1.96*pred.data2$pred.logit.se
pred.data2$logit.U = pred.data2$pred.logit + 1.96*pred.data2$pred.logit.se  
pred.data2$L = exp(pred.data2$logit.L)/(1+exp(pred.data2$logit.L))      
pred.data2$U = exp(pred.data2$logit.U)/(1+exp(pred.data2$logit.U))
pred.data2$U[is.na(pred.data2$U)] =1 
pred.data2$logit.U[pred.data2$logit.U>40] =40
pred.data2$logit.L[pred.data2$logit.L< (-40)] =-40

ggplot(pred.data2,aes(x=Year,y=pred))+  
  facet_wrap(~Age, ncol = 2,scales="free_y") + 
  geom_errorbar(aes(ymin=L,ymax=U),width=0.2,color='blue')+    
  geom_line()+  
  xlab("Year")+ylab("Proportion Mature, with CI") +
  theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
        legend.position="top",legend.key.width=unit(2, 'cm'),
        legend.text = element_text(size=8)) +
  theme_bw()

# ggplot(pred.data2,aes(x=Year,y=pred.logit))+  
#   facet_wrap(~Age, ncol = 2,scales="free_y") +  
#   geom_errorbar(aes(ymin=logit.L,ymax=logit.U),width=0.2,color='grey')+  
#   geom_line()+  
#   xlab("Year")+ylab("Logit Mature, with CI") +
#   theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
#         legend.position="top",legend.key.width=unit(2, 'cm'),
#         legend.text = element_text(size=8)) +
#   theme_bw()

mo[[2]] <- acast(pred.data2, fyear ~ fage, value.var = "pred")


#########################
# Maturity ogive - Viking
#########################

rm(list=ls()[!ls() %in% c("mat_data","ages","years","mo")])

dat = subset(mat_data,SubArea=="V")

dat$wt = dat$wi*dat$wtC
dat$y = dat$NoMature
dat$n = dat$NoAtAlk
dat$fyear = factor(dat$Year)
dat$Year = as.numeric(as.character(dat$fyear))  
dat$fage=factor(dat$Age)
dat$Age = dat$Age

## M4 has a fixed age effects and iid random effects for year and age-year devs. Not really a model but a data summary

m4 <- glmmTMB(cbind(y,n-y) ~ fage - 1+ (1|fyear) + (1|(fage:fyear)),
              data  = dat, weight=wt,family = "binomial")

pred.data2=expand.grid(Age=ages,Year=years)  
pred.data2$fage = factor(pred.data2$Age)        
pred.data2$fyear = factor(pred.data2$Year)       
pred.data2$wt = 1

tmp = predict(m4,type="response",newdata=pred.data2,se.fit=TRUE,allow.new.levels=TRUE)
pred.data2$pred=tmp$fit
pred.data2$pred.se = tmp$se.fit  
tmp = predict(m4,type="link",newdata=pred.data2,se.fit=TRUE,allow.new.levels=TRUE)  
pred.data2$pred.logit=tmp$fit
pred.data2$pred.logit.se = tmp$se.fit

pred.data2$logit.L = pred.data2$pred.logit - 1.96*pred.data2$pred.logit.se
pred.data2$logit.U = pred.data2$pred.logit + 1.96*pred.data2$pred.logit.se  
pred.data2$L = exp(pred.data2$logit.L)/(1+exp(pred.data2$logit.L))      
pred.data2$U = exp(pred.data2$logit.U)/(1+exp(pred.data2$logit.U))
pred.data2$U[is.na(pred.data2$U)] =1 
pred.data2$logit.U[pred.data2$logit.U>40] =40
pred.data2$logit.L[pred.data2$logit.L< (-40)] =-40

ggplot(pred.data2,aes(x=Year,y=pred))+  
  facet_wrap(~Age, ncol = 2,scales="free_y") + 
  geom_errorbar(aes(ymin=L,ymax=U),width=0.2,color='red')+    
  geom_line()+  
  xlab("Year")+ylab("Proportion Mature, with CI") +
  theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
        legend.position="top",legend.key.width=unit(2, 'cm'),
        legend.text = element_text(size=8)) +
  theme_bw()

# ggplot(pred.data2,aes(x=Year,y=pred.logit))+  
#   facet_wrap(~Age, ncol = 2,scales="free_y") +  
#   geom_errorbar(aes(ymin=logit.L,ymax=logit.U),width=0.2,color='grey')+  
#   geom_line()+  
#   xlab("Year")+ylab("Logit Mature, with CI") +
#   theme(strip.text.x = element_text(margin = margin(0.0,0,0.0,0, "cm")),
#         legend.position="top",legend.key.width=unit(2, 'cm'),
#         legend.text = element_text(size=8)) +
#   theme_bw()

mo[[3]] <- acast(pred.data2, fyear ~ fage, value.var = "pred")

names(mo) <- c("south","northwest","viking")
save(mo, file="model/maturity.Rdata")
