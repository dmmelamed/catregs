###
# Updated 6/17/2024 following release of catregs on CRAN
###

rm(list=ls())
require(tidyverse)
# install.packages("catregs")
require(catregs)
data(essUK)
X <- essUK

table(X$walk.alone.dark)
class(X$walk.alone.dark)
table(X$employment)

X <- mutate(X, emp1=ifelse(employment=="Employee",1,0),emp2=ifelse(employment=="Self-employed",1,0),emp3=ifelse(employment=="Unemployed",1,0))
X <- mutate(X,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0))

dat <- X %>% drop_na(walk.alone.dark, minority, female, age, emp1, emp2, religious)

require(MASS)
m1<-polr(walk.alone.dark ~ religious + minority  + female + age + emp1 + emp2,data=dat, Hess=TRUE)
summary(m1)
list.coef(m1)

# Wald test on religious
list.coef(m1)[1,4]^2 # Wald-statistic (z-value squared)
pchisq(list.coef(m1)[1,4]^2,1,lower.tail=FALSE) # p-value on the Wald-test

# LR test on religious
m0<-polr(walk.alone.dark ~  minority  + female + age + emp1 + emp2,data=dat, Hess=TRUE)
anova(m0,m1)

# LR test for multiple coefficients (employment)
m0<-polr(walk.alone.dark ~  religious + minority  + female + age, data=dat, Hess=TRUE)
anova(m1,m0)

# Figure 7.2
list.coef(m1)
m1out <- list.coef(m1)
ggplot(m1out,aes(y=variables,x=b,xmin=ll,xmax=ul)) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=0) +
  labs(x="Log-Odds Coefficient",y="")



# Add reorder option and change variable names; Figure 7.3
ggplot(m1out,aes(y=reorder(variables,b),x=b,xmin=ll,xmax=ul)) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=0) +
  labs(x="Log-Odds Coefficient",y="") +
  scale_y_discrete(labels=c("female"="Female",
                            "emp1"="Traditional Employment",
                            "age"="Age",
                            "religious"="Religious",
                            "minority"="Ethnic Minority",
                            "emp2"="Self-Employed"))


# Odds ratio; Figure 7.4
ggplot(m1out,aes(y=variables,x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b)) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) +
  labs(x="Odds Ratio",y="")
# Add reorder option and change variable names
ggplot(m1out,aes(y=reorder(variables,exp.b),x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b)) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) +
  labs(x="Odds Ratio",y="") + scale_y_discrete(labels=c("female"="Female",
                                                        "emp1"="Traditional Employment",
                                                        "age"="Age",
                                                        "emp2"="Part-Time",
                                                        "religious"="Religious",
                                                        "minority"="Minoritized"))


design <- margins.des(m1,ivs=expand.grid(religious=0:10))
design
pdat <- margins.dat(m1,design)
pdat
ggplot(pdat,aes(x=religious,y=prob,ymin=ll,ymax=ul,group=walk.alone.dark,linetype=walk.alone.dark,color=walk.alone.dark)) +
  theme_bw() + geom_line() + geom_point() + geom_ribbon(alpha=.1) + labs(x="Religious",y="Predicted Probability",color="",linetype="")
require(NatParksPalettes)
?natparks.pals
# Figure 7.5
ggplot(pdat,aes(x=religious,y=prob,ymin=ll,ymax=ul,group=walk.alone.dark,linetype=walk.alone.dark,color=walk.alone.dark)) +
  theme_bw() + geom_line() + geom_point() + geom_ribbon(alpha=.1) + labs(x="Religious",y="Predicted Probability",color="",linetype="") +
  scale_color_manual(values=natparks.pals("Glacier")) +
  theme(legend.position="bottom")


design <- margins.des(m1,ivs=expand.grid(religious=0:10))
pdat <- margins.dat(m1,design)
pdat


# Shaded graph
ggplot(pdat,aes(x=religious,y=prob,fill=walk.alone.dark)) + theme_bw() +
  geom_col(width=.95) +labs(x="Religion",y="Probability",fill="") +
  scale_fill_manual(values=c("grey0","grey30","grey60","grey90")) +
  theme(legend.position="bottom")
ggplot(pdat,aes(x=religious,y=prob,fill=walk.alone.dark)) + theme_bw() +
  geom_col(width=.95) +labs(x="Religion",y="Probability",fill="") +
  scale_fill_manual(values=natparks.pals("Arches",direction=-1)) +
  theme(legend.position="bottom")
# Figure 7.6
ggplot(pdat,aes(x=religious,y=prob,fill=walk.alone.dark)) + theme_bw() +
  geom_col(width=.95) +labs(x="Religion",y="Probability",fill="") +
  scale_fill_manual(values=natparks.pals("Yellowstone")) +
  theme(legend.position="bottom")



# Use cumulative probabilities
pdat <- margins.dat(m1,design,cumulate="yes")
pdat$cut <- factor(pdat$cut,levels=c("Very unsafe|Unsafe",
                                     "Unsafe|Safe",
                                     "Safe|Very safe"))
ggplot(pdat,aes(x=religious,y=cumprob,group=cut,color=cut,ymin=ll,ymax=ul)) +
  theme_bw() + geom_line() + geom_point() + geom_ribbon(alpha=.1) +
  labs(x="Religious",y="Cumulative Probability") +
  theme(legend.position="bottom") +
  scale_color_manual(values=natparks.pals("Banff",direction=-1))

# Figure 7.7
design <- margins.des(m1,ivs=expand.grid(female=c(0,1)))
pdat <- margins.dat(m1,design)
pdat
pdat <- mutate(pdat,xaxs=c(0,.15,.3,.45,.05,.2,.35,.5),
               female=rep(c("Male","Female"),each=4))
ggplot(pdat,aes(x=xaxs,y=prob,ymin=ll,ymax=ul,fill=female)) +
  theme_bw() + geom_col() + labs(x="",y="Predicted Probabilities",fill="") +
  theme(legend.position="bottom") + geom_errorbar(width=.02) +
  scale_x_continuous(breaks=c(.025,.175,.325,.475),labels=c("Very Unsafe","Unsafe","Safe","Very Safe")) +
  scale_fill_manual(values=c("grey65","grey35"))

ggplot(pdat,aes(x=female,y=prob,fill=walk.alone.dark)) + geom_col() +
  labs(x="",y="Probabilities",fill="") + theme_bw() +
  scale_fill_manual(values=c("grey0","grey30","grey60","grey90"))

ggplot(pdat,aes(x=female,y=prob,fill=walk.alone.dark)) + geom_col() +
  labs(x="",y="Probabilities",fill="") + theme_bw() +
  scale_fill_manual(values=natparks.pals("Arches", direction=-1))
# Figure 7.8
ggplot(pdat,aes(x=female,y=prob,fill=walk.alone.dark)) + geom_col() +
  labs(x="",y="Probabilities",fill="") + theme_bw() +
  scale_fill_manual(values=natparks.pals("Cuyahoga"))




# Use cumulative probabilities; Figure 7.9
pdat <- margins.dat(m1,design,cumulate="yes")
pdat <- mutate(pdat,xaxs=c(0,.15,.3,.05,.2,.35),
               female=rep(c("Male","Female"),each=3))
ggplot(pdat,aes(x=xaxs,y=cumprob,ymin=ll,ymax=ul,fill=female)) +
  theme_bw() + geom_col() + labs(x="",y="Predicted Probabilities",fill="") +
  theme(legend.position="bottom") + geom_errorbar(width=.02) +
  scale_x_continuous(breaks=c(.025,.175,.325),labels=c("Very Unsafe|Unsafe","Unsafe|Safe","Safe|Very Safe")) +
  scale_fill_manual(values=c("grey65","grey35"))



# Changes in predicted probabilities
# AMEs
require(marginaleffects)
summary(slopes(m1,variables="female"))
summary(slopes(m1,variables="female",newdata=datagrid(religious=1))) # conditional AME
summary(slopes(m1,variables="female",newdata=datagrid(religious=6))) # conditional AME

mef <- data.frame(slopes(m1,variables="female"))
mef <- mef %>% group_by(group) %>% summarize(estimate=mean(estimate),conf.low=mean(conf.low),conf.high=mean(conf.high))
# Figure 7.10
ggplot(mef,aes(x=group,y=estimate,ymin=conf.low,ymax=conf.high)) + theme_bw() +
  geom_pointrange() + geom_hline(yintercept=0) + labs(x="",y="AME of Female")


mef <- data.frame(slopes(m1,variables="female",newdata=datagrid(religious=1)))
mef2 <- data.frame(slopes(m1,variables="female",newdata=datagrid(religious=6)))
mef<-rbind(mef,mef2)
mef <-mutate(mef,Religion=rep(c("Low","High"),each=4),
             xaxs=c(0,.15,.3,.45,.05,.2,.35,.5))
# Figure 7.11
ggplot(mef,aes(x=xaxs,y=estimate,ymin=conf.low,ymax=conf.high,color=Religion)) + theme_bw() +
  geom_pointrange() + geom_hline(yintercept=0) + labs(x="",y="AME of Female") +
  scale_color_manual(values=c("grey0","grey60")) +
  theme(legend.position="bottom") +
  scale_x_continuous(breaks=c(.025,.175,.325,.475),labels=c("Very unsafe","Unsafe","Safe","Very Safe"))

require(brant)
brant(m1) # overall model violates the parallel lines assumptions


# Partial proportional odds model
dat <- mutate(dat,walk.alone.dark=factor(walk.alone.dark,ordered=TRUE))
require(VGAM)
m2 <- vglm(walk.alone.dark ~ religious + minority  + female + age + emp1 + emp2,
           cumulative(parallel = FALSE ~ age + female, reverse = FALSE),
           data = dat)
summary(m2)
AIC(m1);AIC(m2)
BIC(m1);BIC(m2)
list.coef(m2)
design <- margins.des(m1,ivs=expand.grid(age=seq(18,88,5))) # use the original model to develop the design matrix
design

# Figure 7.13
pdat2 <- margins.dat(m1,design)
pdat2
p1 <- ggplot(pdat2,aes(x=age,y=prob,ymin=ll,ymax=ul,group=walk.alone.dark,linetype=walk.alone.dark,color=walk.alone.dark)) +
  theme_bw() + geom_line() + geom_point() + geom_ribbon(alpha=.1) +
  scale_color_manual(values=natparks.pals("DeathValley")) + labs(x="Age",y="Probability",linetype="",color="")









