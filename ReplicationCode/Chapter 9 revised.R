###
# Updated 5/24/2024. Some functions changed slightly.
### 

rm(list=ls())
require(tidyverse)
require(catregs)
require(marginaleffects)

#Count models
# Poisson
# Interpreting output
# Predicted Counts
# Predicted probabilities
# Residual diagnostics
# Then Negative Binomial
# Zero Inflated Poisson and Negative Binomial
# Count Fit
# Zero Truncated
# Hurdle
data(essUK)

X <- mutate(essUK,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0),
            divorced=ifelse(marital=="Divorced",1,0),married=ifelse(marital=="Married",1,0),widow=ifelse(marital=="Widow",1,0))

m1 <- glm(num.children ~ religious + minority  + female + age + education + divorced + married + widow ,data=X,family="poisson")
summary(m1)
list.coef(m1)
# Figure 9.1
ggplot(X,aes(x=num.children)) + geom_bar() + theme_bw() +
  labs(x="Number of Children",y="Count") + scale_x_continuous(breaks=0:10)




# Coefficient plots; Figure 9.2
pdat <- list.coef(m1)
ggplot(pdat[2:nrow(pdat),],aes(y=reorder(variables,exp.b),x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b)) +
  theme_bw() + geom_pointrange() + labs(x="Incidence Rate Ratio",y="") + geom_vline(xintercept=1)
ggplot(pdat[2:nrow(pdat),],aes(x=b,xmin=ll,xmax=ul,y=variables)) +
  theme_bw() + geom_pointrange() + labs(x="Logged-Coefficient",y="") + geom_vline(xintercept=0)


# Predicted counts; Figure 9.3
design <- margins.des(m1,ivs=expand.grid(education=10:18))
pdat <- margins.dat(m1,design)
ggplot(pdat,aes(x=education,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_point() + geom_line() + geom_ribbon(alpha=.1) +
  labs(x="Education",y="Predicted Count of Children")

first.diff.fitted(m1,pdat,compare=c(1,9))

d1 <- margins.des(m1,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m1,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m1,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m1,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat<- margins.dat(m1,design)
pdat <- mutate(pdat,marital=c("Single","Divorced","Married","Widow"))
# Figure 9.4
ggplot(pdat,aes(x=marital,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_pointrange() + labs(x="",y="Predicted Count of Children")



pdat<- margins.dat(m1,design)
first.diff.fitted(m1,pdat,compare=c(1,2))
first.diff.fitted(m1,pdat,compare=c(1,3))
first.diff.fitted(m1,pdat,compare=c(1,4))

summary(slopes(m1))



# Predicted probabilities; Figure 9.5
table(X$num.children)
pdat <- data.frame(Counts=rep(0:10,2),vals=c(dpois(0:10,mean(m1$fitted.values)),table(X$num.children)/length(X$num.children)),
                   type=rep(c("model","observed"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point() + labs(y="Predicted Probability",linetype="") +
  theme(legend.position = "bottom")


# How do predicted probabilities vary by education?
design <- margins.des(m1,ivs=expand.grid(education=c(10,18)))
pdat <- margins.dat(m1,design)
pdat <- data.frame(Counts=rep(0:10,2),vals=c(dpois(0:10,pdat$fitted[1]),dpois(0:10,pdat$fitted[2])),
                   type=rep(c("10","18"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point()
# How do predicted probabilities vary by marital status?
d1 <- margins.des(m1,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m1,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m1,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m1,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat<- margins.dat(m1,design)
pdat <- data.frame(Counts=rep(0:10,2),
                   vals=c(dpois(0:10,pdat$fitted[1]),dpois(0:10,pdat$fitted[2]),
                          dpois(0:10,pdat$fitted[3]),dpois(0:10,pdat$fitted[4])),
                   type=rep(c("Single","Divorced","Married","Widow"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type,color=type)) +
  theme_bw() + geom_line() + geom_point()


# Diagnostics
# According to Hilbe, the best approach is to plot h-values (hat values) against
# standardized Pearson residuals. Cases high on both are problematic
m1 <- glm(num.children ~ religious + minority  + female + age + education + divorced + married + widow ,data=X,family="poisson")
dia <- diagn(m1)

ggplot() +
  geom_point(data=dia,aes(x=h,y=stdpres),size=.2) + theme_bw() +
  geom_text(data=dia,aes(x=h,y=stdpres,label=obs),position=position_nudge(y=.3)) +
  xlab("Hat Value") + ylab("Standardized Pearson Residual")

# Robustness check
X2 <- data.frame(m1$model,h=dia$h,stdpres=dia$stdpres)
dim(X2)
X2 <- X2 %>% filter(h<.015 | stdpres<2.5)
dim(X2)
m1.check <- glm(num.children ~ religious + minority  + female + age + education + divorced + married + widow ,data=X2,family="poisson")
cbind(coef(m1),coef(m1.check)) # Compare coefficients
cbind(sqrt(diag(vcov(m1))),sqrt(diag(vcov(m1.check)))) # Compare standard errors
cbind(coef(m1)/sqrt(diag(vcov(m1))),coef(m1.check)/sqrt(diag(vcov(m1.check)))) # Compare t-statistics

# Negative binomial
require(MASS)
m2 <- glm.nb(num.children ~ religious + minority  + female + age + education + divorced + married + widow ,data=X)
summary(m2)
lr.test(m1,m2)

# Coefficient plots; Figure 9.6
pdat <- list.coef(m2)
ggplot(pdat[2:nrow(pdat),],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=reorder(variables,exp.b))) +
  theme_bw() + geom_pointrange() + labs(x="Incidence Rate Ratio",y="") + geom_vline(xintercept=1)



ggplot(pdat[2:nrow(pdat),],aes(x=b,xmin=ll,xmax=ul,y=variables)) +
  theme_bw() + geom_pointrange() + labs(x="Logged-Coefficient",y="") + geom_vline(xintercept=0)

# Predicted counts
design <- margins.des(m2,ivs=expand.grid(education=10:18))
pdat <- margins.dat(m2,design)
ggplot(pdat,aes(x=education,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_point() + geom_line() + geom_ribbon(alpha=.1) +
  labs(x="Education",y="Predicted # of Children")

# Figure 9.7
d1 <- margins.des(m2,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m2,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m2,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m2,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat<- margins.dat(m2,design)
pdat <- mutate(pdat,marital=c("Single","Divorced","Married","Widow"))
ggplot(pdat,aes(x=marital,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_pointrange() + labs(x="",y="Predicted Count of Children")


pdat<- margins.dat(m2,design)
first.diff.fitted(m2,pdat,compare=c(1,2))
first.diff.fitted(m2,pdat,compare=c(1,3))
first.diff.fitted(m2,pdat,compare=c(1,4))

summary(marginaleffects(m2))
# Predicted probabilities
pdat <- data.frame(Counts=rep(0:10,2),vals=c(dnbinom(0:10,m2$theta,mu=mean(m2$fitted.values)),table(X$num.children)/length(X$num.children)),
                   type=rep(c("model","observed"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point() + labs(y="Predicted Probability") +
  scale_x_continuous(breaks=0:10) + theme(legend.position="bottom")


# How do predicted probabilities vary by education?
design <- margins.des(m2,ivs=expand.grid(education=c(10,18)))
pdat <- margins.dat(m2,design)
pdat <- data.frame(Counts=rep(0:10,2),vals=c(dnbinom(0:10,m2$theta,mu=pdat$fitted[1]),
                   dnbinom(0:10,m2$theta,mu=pdat$fitted[2])),
                   type=rep(c("10","18"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point()
# How do predicted probabilities vary by marital status?
d1 <- margins.des(m2,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m2,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m2,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m2,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat<- margins.dat(m2,design)
pdat <- data.frame(Counts=rep(0:10,2),
                   vals=c(dnbinom(0:10,m2$theta,mu=pdat$fitted[1]),
                          dnbinom(0:10,m2$theta,mu=pdat$fitted[2]),
                          dnbinom(0:10,m2$theta,mu=pdat$fitted[3]),
                          dnbinom(0:10,m2$theta,mu=pdat$fitted[4])),
                   type=rep(c("Single","Divorced","Married","Widow"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type,color=type)) +
  theme_bw() + geom_line() + geom_point()


# Diagnostics
m2 <- glm.nb(num.children ~ religious + minority  + female + age + education + divorced + married + widow ,data=X)
dia <- diagn(m2)

ggplot() +
  geom_point(data=dia,aes(x=h,y=stdpres),size=.2) + theme_bw() +
  geom_text(data=dia,aes(x=h,y=stdpres,label=obs),position=position_nudge(y=.3)) +
  xlab("Hat Value") + ylab("Standardized Pearson Residual")

# Zero inflated models
require(pscl)
m3 <- zeroinfl(num.children ~ religious + minority  + female + age + education + divorced + married + widow | religious + minority  + female + age + education + divorced + married + widow, dist = "poisson", data = X)
summary(m3)
list.coef(m3)
m4 <- zeroinfl(num.children ~ religious + minority  + female + age + education + divorced + married + widow | religious + minority  + female + age + education + divorced + married + widow, dist = "negbin", data = X)
summary(m4)
list.coef(m4)

# Coefficient plots; Figure 9.9
pdat <- list.coef(m4)
ggplot(pdat[2:9,],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=reorder(variables,b))) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) +
  labs(x="Incidence Rate Ratio",y="")
ggplot(pdat[11:17,],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=reorder(variables,b))) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) +
  labs(x="Odds Ratio",y="")
p1<- ggplot(pdat[2:9,],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=reorder(variables,b))) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) +
  labs(x="Incidence Rate Ratio",y="")
p2<- ggplot(pdat[11:17,],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=reorder(variables,b))) +
  theme_bw() + geom_pointrange() + geom_vline(xintercept=1) +
  labs(x="Odds Ratio",y="")
require(ggpubr)
ggarrange(p1,p2,labels=c("A","B"))




# Predicted counts
design <- margins.des(m4,ivs=expand.grid(education=10:18))
pdat <- margins.dat(m4,design) # does not work, need to specify data
pdat <- margins.dat(m4,design,pscl.data=X)

# Figure 9.10
ggplot(pdat,aes(x=education,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_point() + geom_line() + geom_ribbon(alpha=.1) +
  labs(x="Education",y="Predicted Count of Children")


d1 <- margins.des(m4,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m4,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m4,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m4,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat<- margins.dat(m4,design,pscl.data=X)
pdat <- mutate(pdat,marital=c("Single","Divorced","Married","Widow"))
ggplot(pdat,aes(x=marital,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_pointrange() + labs(x="",y="Predicted # of Children") +
  scale_y_continuous(limits=c(.5,2.5))

# Mems
design <- margins.des(m4,expand.grid(education=c(10,18)))
pdat <- margins.dat(m4,design,pscl.data=X)
compare.margins(margins=pdat$fitted,margins.ses=pdat$se)

# Mems
design <- margins.des(m4,expand.grid(female=c(1,0),education=c(10,18)))
pdat <- margins.dat(m4,design,pscl.data=X)
pdat
compare.margins(margins=pdat$fitted[1:2],margins.ses=pdat$se[1:2])
compare.margins(margins=pdat$fitted[3:4],margins.ses=pdat$se[3:4])

# Ames
ma1<-slopes(m4,variables="female",newdata=datagrid(education=10)) # conditional AME
ma2<-slopes(m4,variables="female",newdata=datagrid(education=18)) # conditional AME
ma<-rbind(ma1,ma2)
ma
compare.margins(margins=ma$estimate,margins.ses=ma$std.error)

# Predicted probabilities; Figure 9.11
apply(m4$model,2,FUN="mean")
design <- margins.des(m4,ivs=expand.grid(female=.548))

predict(m4,newdata=design,type="prob") 

pdat <- data.frame(Counts=rep(0:10,2),vals=c(predict(m4,newdata=design,type="prob"),table(X$num.children)/length(X$num.children)),
                   type=rep(c("model","observed"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point() + scale_x_continuous(breaks=0:10) +
  theme(legend.position="bottom") + labs(linetype="")




# How do predicted probabilities vary by education?
design <- margins.des(m4,ivs=expand.grid(education=c(10,18)))
pdat <- data.frame(Counts=rep(0:10,2),
                   vals=c(predict(m4,newdata=design[1,],type="prob"),
                          predict(m4,newdata=design[2,],type="prob")),
                   type=rep(c("10","18"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point()
# How do predicted probabilities vary by marital status?
d1 <- margins.des(m4,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m4,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m4,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m4,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat <- data.frame(Counts=rep(0:10,2),
                   vals=c(predict(m4,newdata=design[1,],type="prob"),
                          predict(m4,newdata=design[2,],type="prob"),
                          predict(m4,newdata=design[3,],type="prob"),
                          predict(m4,newdata=design[4,],type="prob")),
                   type=rep(c("Single","Divorced","Married","Widow"),each=11))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type,color=type)) +
  theme_bw() + geom_line() + geom_point()


# Diagnostics
dia <- diagn(m4)
ggplot() +
  geom_text(data=dia,aes(x=obs,y=pearsonres,label=obs),position=position_nudge(y=.3)) +
  xlab("Case") + ylab("Pearson Residual") + theme_bw()
# Check robustness to outliers.

# Supplemental model, zero-inflated equation is reduced here.
m5 <- zeroinfl(num.children ~ religious + minority  + female + age + education + divorced + married + widow | female + age + divorced + married + widow, dist = "poisson", data = X)
summary(m5)

# Compare Count models
cf <- count.fit(m1,y.range=0:10)
names(cf)
cf$ic
cf$models
cf$models.pic # Figure 9.12

cf$pic # Figure 9.13


# Zero-truncated counts
X <- mutate(X,ztkids=num.children)
X$ztkids[which(X$ztkids==0)]<-NA
table(X$ztkids)

require(countreg)
m6 <- zerotrunc(ztkids ~ religious + minority  + female + age + education + divorced + married + widow, dist = "poisson", data = X)
summary(m6)
list.coef(m6)
m7 <- zerotrunc(ztkids ~ religious + minority  + female + age + education + divorced + married + widow, dist = "negbin", data = X)
summary(m7)
list.coef(m7)

# Coefficient plots; Figure 9.14
pdat <- list.coef(m7)
ggplot(pdat[2:nrow(pdat),],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=reorder(variables,b))) +
  theme_bw() + geom_pointrange() + labs(x="Incidence Rate Ratios",y="") + geom_vline(xintercept=1)



d1 <- margins.des(m7,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m7,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m7,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m7,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat <- margins.dat(m7,design)
pdat <- mutate(pdat,type=c("Single","Divorced","Married","Widow"))
# Fig 9.15
ggplot(pdat,aes(x=type,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_pointrange() + labs(x="",y="Predicted Count of Children") +
  scale_y_continuous(limits=c(2,3))





# Predicted probabilities
design <- margins.des(m7,ivs=expand.grid(female=.548))
predict(m7,newdata=design,type="prob")
pdat <- data.frame(Counts=rep(1:10,2),vals=c(predict(m7,newdata=design,type="prob"),table(X$ztkids)/length(X$ztkids)),
                   type=rep(c("model","observed"),each=10))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point()

# How do predicted probabilities vary by education?
design <- margins.des(m7,ivs=expand.grid(education=c(10,18)))
pdat <- data.frame(Counts=rep(1:10,2),
                   vals=c(predict(m7,newdata=design[1,],type="prob"),
                          predict(m7,newdata=design[2,],type="prob")),
                   type=rep(c("10","18"),each=10))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type)) +
  theme_bw() + geom_line() + geom_point()
# How do predicted probabilities vary by marital status?
d1 <- margins.des(m7,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m7,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m7,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m7,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat <- data.frame(Counts=rep(1:10,2),
                   vals=c(predict(m7,newdata=design[1,],type="prob"),
                          predict(m7,newdata=design[2,],type="prob"),
                          predict(m7,newdata=design[3,],type="prob"),
                          predict(m7,newdata=design[4,],type="prob")),
                   type=rep(c("Single","Divorced","Married","Widow"),each=10))
ggplot(pdat,aes(x=Counts,y=vals,group=type,linetype=type,color=type)) +
  theme_bw() + geom_line() + geom_point()

# Diagnostics
dia <- diagn(m7)
ggplot() +
  geom_text(data=dia,aes(x=obs,y=pearsonres,label=obs),position=position_nudge(y=.3)) +
  xlab("Case") + ylab("Pearson Residual")
# Check robustness to outliers.


# Hurdle Models
require(pscl)
m8 <- hurdle(num.children ~ religious + minority  + female + age + education + divorced + married + widow , dist = "poisson", data = X)
summary(m8)
list.coef(m8)

m9 <- hurdle(num.children ~ religious + minority  + female + age + education + divorced + married + widow , dist = "negbin", data = X)
summary(m9)
list.coef(m9)

# Coefficient plots; Figure 9.16
pdat <- list.coef(m8)
ggplot(pdat[-c(1,10),],aes(x=exp.b,xmin=ll.exp.b,xmax=ul.exp.b,y=variables)) +
  theme_bw() + geom_pointrange() + labs(x="Exp(b)",y="") + geom_vline(xintercept=1)

ggplot(pdat[-c(1,10),],aes(x=b,xmin=ll,xmax=ul,y=variables)) +
  theme_bw() + geom_pointrange() + labs(x="Logged-Coefficient",y="") + geom_vline(xintercept=0)

# Predicted counts; Figure 9.17
design <- margins.des(m8,ivs=expand.grid(education=10:18))
pdat <- margins.dat(m8,design)
ggplot(pdat,aes(x=education,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_point() + geom_line() + geom_ribbon(alpha=.2) +
  labs(x="Education",y="Predicted # of Children")

d1 <- margins.des(m8,ivs=expand.grid(divorced=0,married=0,widow=0))
d2 <- margins.des(m8,ivs=expand.grid(divorced=1,married=0,widow=0))
d3 <- margins.des(m8,ivs=expand.grid(divorced=0,married=1,widow=0))
d4 <- margins.des(m8,ivs=expand.grid(divorced=0,married=0,widow=1))
design <- rbind(d1,d2,d3,d4)
design
pdat<- margins.dat(m8,design)
pdat <- mutate(pdat,marital=c("Single","Divorced","Married","Widow"))
ggplot(pdat,aes(x=marital,y=fitted,ymin=ll,ymax=ul)) +
  theme_bw() + geom_pointrange() + labs(x="",y="Predicted # of Children")

# Predicted probabilities; Figure 9.18
table(X$num.children)
design <- margins.des(m8,ivs=expand.grid(female=.548))

