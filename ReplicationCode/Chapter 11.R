###
# Updated 6/17/2024 following release of catregs on CRAN
###

rm(list=ls())
require(tidyverse)
# install.packages("catregs")
require(catregs)
data(essUK)
X <- filter(essUK,country=="United Kingdom")
# create variables
X <- mutate(X,safe = ifelse(walk.alone.dark=="Safe" | walk.alone.dark==
                              "Very safe",1,0))
X <- mutate(X, employee=ifelse(employment=="Employee",1,0),self.emp=ifelse(employment=="Self-employed",1,0),emp3=ifelse(employment=="Unemployed",1,0))
X <- mutate(X,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0))
X <- mutate(X,imm.good=ifelse(immigration.good.economy>5,1,0))
X2 <- X %>% drop_na(imm.good , religious ,conservative ,  female ,  minority , age,  employment)
# Reduced model
m1 <- glm(imm.good ~  minority +  religious +   female +   age+  self.emp + employee,data=X2,family=binomial)
summary(m1)
# Full model
m2 <- glm(imm.good ~ minority + religious +  female +   age+  self.emp + employee +  conservative,data=X2,family=binomial)
summary(m2)
require(systemfit)
require(sandwich)
k1 <- khb(m1,m2) # load the function. It's available in the replication folder on Github too
print.khb(k1)




m3 <- lm(conservative ~ minority +  religious +  female +   age+  self.emp + employee,data=X2)
X2 <- mutate(X2,resids=m3$residuals)
# Adjusted reduced model
m1.adj <- glm(imm.good ~  minority +  religious +  female +   age+  self.emp + employee + resids,data=X2,family=binomial)
summary(m1.adj) # Coefficient for minority is the "reduced" coefficient in the KHB output

require(margins)
summary(margins(m1.adj))
summary(margins(m2))

margins.reduced<-summary(margins(m1.adj))
margins.full <- summary(margins(m2))
margins.reduced
margins.full
compare.margins(margins=c(margins.reduced$AME[1],margins.full$AME[1]),
                margins.ses=c(margins.reduced$SE[1],margins.full$SE[1])) # Age

compare.margins(margins=c(margins.reduced$AME[2],margins.full$AME[3]),
                margins.ses=c(margins.reduced$SE[2],margins.full$SE[3])) # Employee

compare.margins(margins=c(margins.reduced$AME[3],margins.full$AME[4]),
                margins.ses=c(margins.reduced$SE[3],margins.full$SE[4])) # Female

compare.margins(margins=c(margins.reduced$AME[4],margins.full$AME[5]),
                margins.ses=c(margins.reduced$SE[4],margins.full$SE[5])) # Minority

compare.margins(margins=c(margins.reduced$AME[5],margins.full$AME[6]),
                margins.ses=c(margins.reduced$SE[5],margins.full$SE[6])) # Religious

compare.margins(margins=c(margins.reduced$AME[7],margins.full$AME[7]),
                margins.ses=c(margins.reduced$SE[7],margins.full$SE[7])) # Self-Employed



# Continuous variable
design <- margins.des(m1.adj,ivs=expand.grid(age=seq(18,80,5)))
mar1 <- margins.dat(m1.adj,design)
design2 <- margins.des(m2,ivs=expand.grid(age=seq(18,80,5)))
mar2 <- margins.dat(m2,design)
dim(mar1)
dim(mar2)
# Compare 18 year old to 78 year olds
fd1 <- first.diff.fitted(m1.adj,design,compare=c(13,1))
fd2 <- first.diff.fitted(m2,design2,compare=c(13,1))
compare.margins(margins=c(fd1$`First Difference`,fd2$`First Difference`),
                margins.ses=c(fd1$`Standard Error`,fd2$`Standard Error`),nsim=1000000)

pdat <- rbind(fd1,fd2)
pdat <- mutate(pdat,mod=c("Reduced: Controlling for Residual of Conservative","Full: Controlling for Conservative"))
pdat$mod <- factor(pdat$mod,levels=c("Reduced: Controlling for Residual of Conservative","Full: Controlling for Conservative"))
# Figure 11.1
ggplot(pdat,aes(x=mod,y=`First Difference`,ymin=ll,ymax=ul)) + theme_bw() +
  geom_point() + geom_errorbar(width=.05) + scale_y_continuous(limits=c(-.25,0)) +
  geom_hline(yintercept=0) + labs(x="",y="Marginal Effect at Means (MEM) of Age")




# Missing Data
rm(list=ls())
require(tidyverse)
require(catregs)
data(essUK)
X <- filter(essUK,country=="United Kingdom")
# create variables
X <- mutate(X,safe = ifelse(walk.alone.dark=="Safe" | walk.alone.dark==
                              "Very safe",1,NA))
X$safe[which(X$walk.alone.dark=="Unsafe" | X$walk.alone.dark=="Very unsafe")] <- 0
X <- mutate(X, employee=ifelse(employment=="Employee",1,0),self.emp=ifelse(employment=="Self-employed",1,0),emp3=ifelse(employment=="Unemployed",1,0))

table(X$emp1)
table(X$emp2)
table(X$emp3)

X <- mutate(X,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0))


m1 <- glm(safe ~ religious + minority  + female + age + self.emp + employee,data=X,family=binomial)
summary(m1)
dim(m1$model)
dim(X)

sum(is.na(X$safe))
variables<-c("Feel Safe at Night","Religiousness","Racially Minoritized",
             "Female","Age","Employee","Self-Employed","Unemployed") #Variable names

means <- c(mean(X$safe,na.rm=TRUE),
           mean(X$religious,na.rm=TRUE),
           mean(X$minority,na.rm=TRUE),
           mean(X$female,na.rm=TRUE),
           mean(X$age,na.rm=TRUE),
           table(X$employment)[1]/sum(table(X$employment)),
           table(X$employment)[2]/sum(table(X$employment)),
           table(X$employment)[3]/sum(table(X$employment))) # Variable means/proportions
sds <- c("",
         round(sd(X$religious,na.rm=TRUE),3),
         "",
         "",
         round(sd(X$age,na.rm=TRUE),3),
         "","","") # Standard deviations for quant variables
ns<-c(sum(1-is.na(X$safe)),
      sum(1-is.na(X$religious)),
      sum(1-is.na(X$minority)),
      sum(1-is.na(X$female)),
      sum(1-is.na(X$age)),
      sum(1-is.na(X$employment)),
      sum(1-is.na(X$employment)),
      sum(1-is.na(X$employment))) # sample size per variable
descriptives<-data.frame(variables,means,sds,ns) # put them together
descriptives
descriptives <- mutate(descriptives,mis=1-ns/2204)
# print as table
require(kableExtra)
descriptives %>%
  kbl(digits=3,col.names=c("Variable","Mean","SD","N","% Missing")) %>%
  kable_classic(full_width = F)

require(mice)
m1 <- glm(safe ~ religious + minority  + female + age + self.emp + employee,data=X,family=binomial)


X2 <- X %>% select(safe , religious , minority  , female , age ,self.emp, employee)
X2$employee <- as.factor(X2$employee)
X2$minority<-as.factor(X2$minority)
X2$female<-as.factor(X2$female)
out<-mice(X2,maxit=0)
names(out)
out$method # Default imputation methods
out$predictorMatrix # Rows represent the variables used to impute missing values
predMat <- out$predictorMatrix
predMat[2,c(3,5)]<-0 # now minorty and age are not used in imputing religious
out<-mice(X2,m=20,seed=1982,method=c("","norm","logreg","logreg","norm","norm","logreg"),printFlag = FALSE)
fits<-with(out,glm(safe ~ religious + minority  + female + age + self.emp + employee,data=X,family=binomial))
est<-pool(fits)
summary(est)

require(margins)
margs<-data.frame(summary(margins(fits$analyses[[1]],variables=c("female"))))
for(i in 2:20){
  margsi<-data.frame(summary(margins(fits$analyses[[i]],variables=c("female"))))
  margs<-rbind(margs,margsi)}
margs %>% summarize(AME=mean(AME),SE=rubins.rule(SE))
dnorm(-.1801204/.01853468)

margs<-data.frame(summary(margins(fits$analyses[[1]],variables=c("age"))))
for(i in 2:20){
  margsi<-data.frame(summary(margins(fits$analyses[[i]],variables=c("age"))))
  margs<-rbind(margs,margsi)}
margs %>% summarize(AME=mean(AME),SE=rubins.rule(SE))
dnorm(-0.001226589/.0005211336)




