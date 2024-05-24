###
# Updated 5/24/2024. Some functions changed slightly.
### 

rm(list=ls())
require(tidyverse)
require(catregs)
data(essUK)

X <- filter(essUK,country=="United Kingdom")
table(X$walk.alone.dark)

X <- mutate(X,safe = ifelse(walk.alone.dark=="Safe" | walk.alone.dark==
                              "Very safe",1,0))

X <- mutate(X,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0),
            divorced=ifelse(marital=="Divorced",1,0),married=ifelse(marital=="Married",1,0),widow=ifelse(marital=="Widow",1,0))
X %>% drop_na(minority) %>% group_by(minority,female) %>% summarise(saf=mean(safe,na.rm=TRUE))
pd<-X %>% drop_na(minority) %>% group_by(minority,female) %>% summarise(saf=mean(safe,na.rm=TRUE))
pd <- mutate(data.frame(pd),cond=c("Male/Non-Minority","Female/Non-Minority",
                                   "Male/Minority","Female/Minority"))
pd$cond<-factor(pd$cond,levels=c("Female/Non-Minority","Male/Non-Minority","Female/Minority","Male/Minority"))
# Figure 6.1
ggplot(pd,aes(x=cond,y=saf)) + theme_bw() + geom_col() +
  labs(x="",y="Proportion Feeling Safe at Night") + scale_y_continuous(limits=c(0,1))

###
# Categorical x categorical
###
m1 <- glm(safe ~ religious + ethnic.minority  + gender + age ,data=X,family="binomial")
summary(m1)
m2 <- glm(safe ~ religious + minority*female + age,data=X,family="binomial")
design <- margins.des(m2,ivs=expand.grid(minority=c(0,1),female=c(0,1)))
pdat <- margins.dat(m2,design)
pdat
pdat <- mutate(pdat,Minority=rep(c("No","Yes"),2),
               Female=rep(c("Male","Female"),each=2),
               xaxs=c(-.05,0,.1,.15))
# Figure 6.2
ggplot(pdat,aes(x=xaxs,y=fitted,ymin=ll,ymax=ul,group=Minority,color=Minority)) +
  geom_pointrange() + theme_bw() + labs(x="",y="Pr(Safe at night)") +
  scale_x_continuous(breaks=c(-.025,.125),labels=c("Male","Female"),limits=c(-.1,.2)) +
  theme(legend.position="bottom") +
  scale_color_manual(values=c("grey0","grey60"))

  
pdat
first.diff.fitted(m2,design,compare=c(3,1,4,2))

first.diff.fitted(m2,design,compare=c(3,1)) # Effect of gender for racial majority members
first.diff.fitted(m2,design,compare=c(4,2)) # Effect of gender for racial minority members
second.diff.fitted(m2,design,compare=c(3,1,4,2)) # The effect of gender is larger for majority members
mem1<-first.diff.fitted(m2,design,compare=c(3,1))
mem2<-first.diff.fitted(m2,design,compare=c(4,2))
compare.margins(margins=c(mem1$`First Difference`,mem2$`First Difference`),margins.ses=c(mem1$`Standard Error`,mem2$`Standard Error`))


require(margins) # install from archive off of CRAN
ma1 <- summary(margins(m2,variables="female",at=list(minority=0)))
ma2 <- summary(margins(m2,variables="female",at=list(minority=1)))
cames <- rbind(ma1,ma2)
cames # These are conditional average marginal effects
compare.margins(margins=cames$AME,margins.ses=cames$SE)

###
# Categorical x continuous
###

head(X)

m1 <- glm(safe ~ religious + minority  + female  + immigration.good.economy + age ,data=X,family="binomial")
summary(m1)
m2 <- glm(safe ~ religious + minority  + female*immigration.good.economy + age ,data=X,family="binomial")
summary(m2)

design <- margins.des(m2,ivs=expand.grid(female=c(0,1),immigration.good.economy=0:10))
pdat <- margins.dat(m2,design)
pdat <- mutate(pdat,Female=rep(c("Male","Female"),11))
# Figure 6.3
ggplot(pdat,aes(x=immigration.good.economy,y=fitted,ymin=ll,ymax=ul,group=Female,color=Female)) +
  theme_bw() + geom_line() + geom_point() + geom_ribbon(alpha=.2) + labs(x="Immigration good for the economy",y="Pr(Safe at night)",color="") +
  theme(legend.position="bottom") + scale_x_continuous(breaks=seq(0,10,2)) +
  scale_color_manual(values=c("grey0","grey60"))


first.diff.fitted(m2,design,compare=c(2,1,22,21)) # Effect of gender for for those with immigration==0
first.diff.fitted(m2,design,compare=c(22,21)) # Effect of gender for for those with immigration==10
second.diff.fitted(m2,design,compare=c(2,1,22,21)) # The effect of gender is larger when immigration==0 vs 10
second.diff.fitted(m2,design,compare=c(10,9,20,19)) # The effect of gender is larger when immigration==4 vs 9

ma1 <- summary(margins(m2,variables="immigration.good.economy",at=list(female=0)))
ma1
ma2 <- summary(margins(m2,variables="immigration.good.economy",at=list(female=1)))
ma2
compare.margins(margins=c(ma2$AME,ma1$AME),margins.ses=c(ma2$SE,ma1$SE))


design <- margins.des(m2,ivs=expand.grid(immigration.good.economy=0:10,female=c(0,1)))
design
fdm<-first.diff.fitted(m2,design,compare=c(11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,1))
fdf<-first.diff.fitted(m2,design,compare=c(22,21,21,20,20,19,19,18,18,17,17,16,16,15,15,14,14,13,13,12))
mean(fdf$`First Difference`)
rubins.rule(fdf$`Standard Error`)
dnorm(mean(fdf$`First Difference`)/rubins.rule(fdf$`Standard Error`))

fdm
mean(fdm$`First Difference`)
rubins.rule(fdm$`Standard Error`)
dnorm(mean(fdm$`First Difference`)/rubins.rule(fdm$`Standard Error`))



sd1<-second.diff.fitted(m2,design,compare=c(22,21,11,10)) # The effect of gender is larger when immigration==0 vs 10
sd2<-second.diff.fitted(m2,design,compare=c(21,20,10,9)) # The effect of gender is larger when immigration==0 vs 10
sd3<-second.diff.fitted(m2,design,compare=c(20,19,9,8)) # The effect of gender is larger when immigration==0 vs 10
sd4<-second.diff.fitted(m2,design,compare=c(19,18,8,7)) # The effect of gender is larger when immigration==0 vs 10
sd5<-second.diff.fitted(m2,design,compare=c(18,17,7,6)) # The effect of gender is larger when immigration==0 vs 10
sd6<-second.diff.fitted(m2,design,compare=c(17,16,6,5)) # The effect of gender is larger when immigration==0 vs 10
sd7<-second.diff.fitted(m2,design,compare=c(16,15,5,4)) # The effect of gender is larger when immigration==0 vs 10
sd8<-second.diff.fitted(m2,design,compare=c(15,14,4,3)) # The effect of gender is larger when immigration==0 vs 10
sd9<-second.diff.fitted(m2,design,compare=c(14,13,3,2)) # The effect of gender is larger when immigration==0 vs 10
sd10<-second.diff.fitted(m2,design,compare=c(13,12,2,1)) # The effect of gender is larger when immigration==0 vs 10
secdiffs<-rbind(sd1,sd2,sd3,sd4,sd5,sd6,sd7,sd8,sd9,sd10)
mean(secdiffs$`Second Difference`)
rubins.rule(secdiffs$`Standard Error`)
dnorm(mean(secdiffs$`Second Difference`)/rubins.rule(secdiffs$`Standard Error`))












###
# Continuous x continuous
###

m1 <- glm(safe ~ religious + minority  + female  + immigration.good.economy + age ,data=X,family="binomial")
summary(m1)
m2 <- glm(safe ~ religious + minority  + female + religious*immigration.good.economy + age ,data=X,family="binomial")
summary(m2)

table(X$religious)
table(X$immigration.good.economy)

design <- margins.des(m2,ivs=expand.grid(immigration.good.economy=c(0,3,6),
                      religious=0:10))
pdat <- margins.dat(m2,design)

# Not Figure 6.4
ggplot(pdat,aes(x=religious,y=fitted,ymin=ll,ymax=ul,group=immigration.good.economy,color=immigration.good.economy)) +
  theme_bw() + geom_point() + geom_line() + scale_color_gradient(low="purple",high="orange") +
  labs(color="Immigration good for the economy",x="Religiousness",y="Pr(Safe at Night)") +
  theme(legend.position="bottom") + scale_x_continuous(breaks=0:10)
design <- margins.des(m2,ivs=expand.grid(immigration.good.economy=0:10,
                                         religious=c(0,3,6)))
pdat <- margins.dat(m2,design)
# Figure 6.4
ggplot(pdat,aes(x=immigration.good.economy,y=fitted,ymin=ll,ymax=ul,group=religious,color=religious)) +
  theme_bw() + geom_point() + geom_line() + scale_color_gradient(low="purple",high="orange") +
  labs(color="Religiousness",x="Immigration good for the economy",y="Pr(Safe at Night)") +
  theme(legend.position="bottom") + scale_x_continuous(breaks=0:10)


design <- margins.des(m2,ivs=expand.grid(immigration.good.economy=0:10,
                                         religious=0))
fd1<-first.diff.fitted(m2,design,compare=c(11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,1))
fd1<-mutate(fd1,relig=0)
for(i in 1:10){
  design <- margins.des(m2,ivs=expand.grid(immigration.good.economy=0:10,
                                           religious=i))
  fd2<-first.diff.fitted(m2,design,compare=c(11,10,10,9,9,8,8,7,7,6,6,5,5,4,4,3,3,2,2,1))
  fd2<-mutate(fd2,relig=i)
  fd1<-rbind(fd1,fd2)}

fds<- fd1 %>% group_by(relig) %>% summarize(fds=mean(`First Difference`),ses=rubins.rule(`Standard Error`))
fds
compare.margins(margins=fds$fds[c(1,11)],margins.ses=fds$ses[c(1,11)])

# Figure 6.5
ggplot(fds,aes(x=relig,y=fds,ymin=fds-1.96*ses,ymax=fds+1.96*ses)) +
  theme_bw() + geom_line() + geom_ribbon(alpha=.2) + geom_hline(yintercept=0,linetype=2) +
  scale_x_continuous(breaks=0:10) + labs(y="cMEM of Immigration on Feeling Safe",x="Religiousness")


# Figure 6.6
ma1 <- summary(margins(m2,variables="immigration.good.economy",at=list(religious=0)))
for(i in 1:10){
  ma2 <- summary(margins(m2,variables="immigration.good.economy",at=list(religious=i)))
  ma1<-rbind(ma1,ma2)}
ma1[,-1]
ggplot(ma1,aes(x=religious,y=AME,ymin=lower,ymax=upper)) +
  theme_bw() + geom_line() + geom_ribbon(alpha=.2) + geom_hline(yintercept=0,linetype=2) +
  scale_x_continuous(breaks=0:10) + labs(y="cAME of Immigration on Feeling Safe",x="Religiousness")

compare.margins(margins=ma1$AME[c(1,11)],margins.ses=ma1$SE[c(1,11)])

# Contour plot/heatmap, Figure 6.7
design <- margins.des(m2,ivs=expand.grid(immigration.good.economy=0:10,
                                         religious=0:10))
pdat <- margins.dat(m2,design)
ggplot(pdat,mapping=aes(x=immigration.good.economy, y=religious, fill=fitted)) +
  geom_tile() + theme_bw() + labs(x="Immigration is good for economy",y="Religiousness",fill="Pr(Safe at Night)") +
  scale_x_continuous(breaks=seq(0,10,2)) +scale_y_continuous(breaks=seq(0,10,2)) + theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange")


###
# Squared terms
###
X <- mutate(X,highinc = ifelse(income.decile>6,1,0),
            age2 = age*age)

X2 <- X %>% drop_na(highinc, religious , minority  , female  ,immigration.good.economy , age)
m2 <- glm(highinc ~ religious + minority  + female  + married + age + I(age^2) ,data=X2,family="binomial")
summary(m2)
m1 <- glm(highinc ~ religious + minority  + female  + married + age ,data=X2,family="binomial")
summary(m1)
# Figure 6.8
design<-margins.des(m1,expand.grid(age=20:80))
mar1<-margins.dat(m1,design)
mar2<-margins.dat(m2,design)
dim(mar1)
mars <- rbind(mar1,mar2)
mars<-mutate(mars,type=rep(c("Linear","Squared"),each=nrow(mar1)))
ggplot(mars,aes(x=age,y=fitted,ymin=ll,ymax=ul,group=type,fill=type)) + theme_bw() +
  geom_line() + geom_ribbon(alpha=.4) + theme(legend.position="bottom") +
  labs(y="Pr(High Income)",x="Age",fill="Effect of age:") +
  scale_fill_manual(values=c("purple","orange"))


design2<-margins.des(m2,expand.grid(age=30:80))
design2$`I(age^2)` <- design2$age^2
mems1<-first.diff.fitted(m2,design2,compare=c(2,1))
for (i in 2:(nrow(design2) -1)){
  mems2<-first.diff.fitted(m2,design2,compare=c(i+1,i))
  mems1<-rbind(mems1,mems2)}
mean(mems1$`First Difference`)
rubins.rule(mems1$`Standard Error`)
dnorm(mean(mems1$`First Difference`)/rubins.rule(mems1$`Standard Error`))

summary(margins(m2,variables="age"))
#summary(margins(m1))


design2<-margins.des(m2,expand.grid(age=20:80))
design2$`I(age^2)` <- design2$age^2
mems1<-first.diff.fitted(m2,design2,compare=c(2,1))
for (i in 2:(nrow(design2) -1)){
  mems2<-first.diff.fitted(m2,design2,compare=c(i+1,i))
  mems1<-rbind(mems1,mems2)}
mems1<-mutate(mems1,age=20:79)
mean(mems1$`First Difference`[1:26])
rubins.rule(mems1$`Standard Error`[1:26])
dnorm(mean(mems1$`First Difference`[1:26])/rubins.rule(mems1$`Standard Error`[1:26]))

mean(mems1$`First Difference`[31:60])
rubins.rule(mems1$`Standard Error`[31:60])
dnorm(mean(mems1$`First Difference`[31:60]),rubins.rule(mems1$`Standard Error`[31:60]))
compare.margins(margins=c(mean(mems1$`First Difference`[1:26]),mean(mems1$`First Difference`[31:60])),
                margins.ses=c(rubins.rule(mems1$`Standard Error`[1:26]),rubins.rule(mems1$`Standard Error`[31:60])))


m3 <- glm(highinc ~ religious + female  + minority + married*age + married*I(age^2) ,data=X,family="binomial")
summary(m3)
# Figure 6.9
design <- margins.des(m3,ivs=expand.grid(age=seq(25,65,5),married=c(1,0)))
design[,ncol(design)] <- design$age^2
pdat <- margins.dat(m3,design)
pdat <- mutate(pdat,Married=rep(c("Married","Not Married"),each=nrow(pdat)/2))
ggplot(pdat,aes(x=age,y=fitted,ymin=ll,ymax=ul,group=Married,fill=Married)) +
  theme_bw() + geom_point() + geom_line() + geom_ribbon(alpha=.4) +
  scale_y_continuous(limits=c(-.1,.8)) +
  labs(x="Age",y="Pr(High Wage)",color="") + scale_x_continuous(limits=c(25,65),breaks=seq(30,80,10)) +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("purple","orange"))


# Figure 6.10
design <- margins.des(m3,ivs=expand.grid(married=c(1,0),age=25:65))
design[,ncol(design)] <- design$age^2 #fix age-squared
head(design)
dim(design)
fds <- first.diff.fitted(m3,design,compare=1:nrow(design))
fds 
fds<-mutate(fds,age=25:65)
ggplot(fds,aes(x=age,y=`First Difference`,ymin=ll,ymax=ul))  + theme_bw() + geom_line() +
  geom_ribbon(alpha=.2) + geom_hline(yintercept=0,linetype=2) +
  labs(x="Age",y="MEM of Married, conditional on Age") + scale_y_continuous(limits=c(-.65,.14))

p1<-ggplot(fds,aes(x=age,y=`First Difference`,ymin=ll,ymax=ul))  + theme_bw() + geom_line() +
  geom_ribbon(alpha=.2) + geom_hline(yintercept=0,linetype=2) +
  labs(x="Age",y="MEM of Married, conditional on Age") + scale_y_continuous(limits=c(-.65,.14))



ma1<-summary(margins(m3,variables=c("married"),at=list(age=25)))
for(i in 26:65){
  ma2<-summary(margins(m3,variables=c("married"),at=list(age=i)))
  ma1<-rbind(ma1,ma2)}

ggplot(ma1,aes(x=age,y=AME,ymin=lower,ymax=upper)) + theme_bw() + geom_line() +
  geom_ribbon(alpha=.2) + geom_hline(yintercept=0,linetype=2) +
  labs(x="Age",y="AME of Married, conditional on Age") + scale_y_continuous(limits=c(-.65,.14))
p2<-ggplot(ma1,aes(x=age,y=AME,ymin=lower,ymax=upper)) + theme_bw() + geom_line() +
  geom_ribbon(alpha=.2) + geom_hline(yintercept=0,linetype=2) +
  labs(x="Age",y="AME of Married, conditional on Age") + scale_y_continuous(limits=c(-.65,.14))
require(ggpubr)
ggarrange(p1,p2,labels=c("A","B"))

# Compare 40 to 65
ma1[c(16,41),]
fds[c(16,41),]
compare.margins(margins=ma1$AME[c(16,41)],margins.ses=ma1$SE[c(16,41)])
compare.margins(margins=fds$`First Difference`[c(16,41)],margins.ses=fds$`Standard Error`[c(16,41)])



ma1<-summary(margins(m3,variables=c("age"),at=list(married=0)))
ma2<-summary(margins(m3,variables=c("age"),at=list(married=1)))
ma1;ma2
compare.margins(margins=c(ma1$AME,ma2$AME),margins.ses=c(ma1$SE,ma2$SE))



