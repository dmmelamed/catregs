
rm(list=ls())
require(tidyverse)
require(catregs)
data(essUK)
head(essUK)

X <- mutate(essUK,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0),
            divorced=ifelse(marital=="Divorced",1,0),married=ifelse(marital=="Married",1,0),widow=ifelse(marital=="Widow",1,0))
X <- mutate(X,selfemp=ifelse(employment=="Self-employed",1,0),
            unemp=ifelse(employment=="Unemployed",1,0))

require(nnet)
m1<-multinom(walk.alone.dark ~ education + religious + minority  + female + age + selfemp + unemp,data=X)
summary(m1)
list.coef(m1)

# Change reference category
X$walk.alone.dark2 <- relevel(as.factor(X$walk.alone.dark), ref = "Very safe")
m1<-multinom(walk.alone.dark2 ~ education + religious + minority  + female + age + selfemp + unemp,data=X)
summary(m1)
list.coef(m1)
pdat<-list.coef(m1)$out
pdat <- pdat[-c(1,9,17),]

#Coefficient plot; Figure 8.1
ggplot(pdat,aes(x=b,xmin=ll,xmax=ul,y=variables)) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=0) + labs(y="")



#OR plot; Figure 8.2
ggplot(pdat,aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=variables)) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) + labs(y="")

# Fig 8.3
require(NatParksPalettes)
design <- margins.des(m1,ivs=expand.grid(education=10:18), data=X)
pdat <- margins.dat(m1,design)
ggplot(pdat,aes(x=education,y=prob,ymin=ll,ymax=ul,group=walk.alone.dark2,linetype=walk.alone.dark2,color=walk.alone.dark2)) +
  theme_bw() + geom_line() + geom_point() + geom_ribbon(alpha=.1) +
  labs(x="Education",y="Predicted Probability",linetype="",color="") +
  scale_color_manual(values=natparks.pals("Glacier")) +
  theme(legend.position="bottom")

first.diff.fitted(m1,design,compare=c(1,9))
first.diff.fitted(m1,design,compare=c(3,7))






# Figure 8.4
design <- margins.des(m1,ivs=expand.grid(female=c(0,1)), data=X)
pdat <- margins.dat(m1,design)
pdat <- mutate(pdat,sex=rep(c("Male","Female"),each=4),
               xaxs=c(0,.15,.3,.45,.05,.2,.35,.5))
ggplot(pdat,aes(x=xaxs,y=prob,ymin=ll,ymax=ul,fill=sex)) +
  theme_bw() + geom_col() +
  scale_fill_manual(values=c("grey65","grey35")) +
  geom_errorbar(width=.02) + theme(legend.position="bottom") +
  labs(x="",y="Predicted Probability",fill="") +
  scale_x_continuous(breaks=c(.025,.175,.325,.525),
                     labels=c("Very Unsafe","Unsafe","Safe","Very Safe"))

first.diff.fitted(m1,design,compare=c(1,2)) # Significance of those differences




# Average Marginal Effects and Conditional Average Marginal Effects
require(marginaleffects)
summary(marginaleffects(m1))


summary(marginaleffects(m1,newdata=datagrid(minority=1)))
summary(marginaleffects(m1,newdata=datagrid(minority=0)))

mef <- data.frame(summary(marginaleffects(m1,variables="female",newdata=datagrid(minority=0))))
mef2 <- data.frame(summary(marginaleffects(m1,variables="female",newdata=datagrid(minority=1))))
mef<-rbind(mef,mef2)
mef <-mutate(mef,Minority=rep(c("No","Yes"),each=4),
             xaxs=c(0,.15,.3,.45,.05,.2,.35,.5))
# Figure 8.5
ggplot(mef,aes(x=xaxs,y=estimate,ymin=conf.low,ymax=conf.high,color=Minority)) + theme_bw() +
  geom_pointrange() + geom_hline(yintercept=0) + labs(x="",y="AME of Female") +
  scale_color_manual(values=c("grey0","grey60")) + 
  theme(legend.position="bottom") +
  scale_x_continuous(breaks=c(.025,.175,.325,.475),labels=c("Very unsafe","Unsafe","Safe","Very Safe"))


# An LR test for combining alternatives
m1<-multinom(walk.alone.dark ~ education + religious + minority  + female + age + selfemp + unemp,data=X)
X <- mutate(X,walk2=recode_factor(walk.alone.dark,"Very unsafe"="Unsafe"))
m2<-multinom(walk2 ~ education + religious + minority  + female + age + selfemp + unemp,data=X)
lr.test(m1,m2)



# IIA test
require(mlogit)
X2<-mlogit.data(X,choice="walk.alone.dark",shape="wide")
ml1 <- mlogit(walk.alone.dark ~ 0|education + religious + minority  + female + age + selfemp + unemp, data=X2, reflevel = "Very safe")
ml2 <- mlogit(walk.alone.dark ~ 0|education + religious + minority  + female + age + selfemp + unemp, data=X2, 
              reflevel = "Very safe",
              alt.subset=c("Safe","Unsafe","Very safe"))
hmftest(ml1,ml2) # Test for Very unsafe
ml3 <- mlogit(walk.alone.dark ~ 0|education + religious + minority  + female + age + selfemp + unemp, data=X2, 
              reflevel = "Very safe",
              alt.subset=c("Safe","Very unsafe","Very safe"))
hmftest(ml1,ml3) # Test for Unsafe
ml4 <- mlogit(walk.alone.dark ~ 0|education + religious + minority  + female + age + selfemp + unemp, data=X2, 
              reflevel = "Very safe",
              alt.subset=c("Unsafe","Very unsafe","Very safe"))
hmftest(ml1,ml4) # Test for Safe
ml5 <- mlogit(walk.alone.dark ~ 0|education + religious + minority  + female + age + selfemp + unemp, data=X2, 
              reflevel = "Safe",
              alt.subset=c("Unsafe","Very unsafe","Safe"))
hmftest(ml1,ml5) # Test for Very safe








