rm(list=ls())
require(tidyverse)
require(catregs)
data(essUK)

X <- filter(essUK,country=="United Kingdom")
table(X$walk.alone.dark)

# create variables
X <- mutate(X,safe = ifelse(walk.alone.dark=="Safe" | walk.alone.dark==
                              "Very safe",1,0))

X <- mutate(X, emp1=ifelse(employment=="Employee",1,0),emp2=ifelse(employment=="Self-employed",1,0),emp3=ifelse(employment=="Unemployed",1,0))

X <- mutate(X,minority=ifelse(ethnic.minority=="Yes",1,0),female=ifelse(gender=="Female",1,0))

dat <- X %>% drop_na(safe, minority, female, age, emp1, emp2, religious)

# basic model
m1 <- glm(safe ~ religious + minority  + female + age + emp1 + emp2,data=dat,family=binomial)
summary(m1) # Table 5.1
m1probit <- glm(safe ~ religious + minority + female + age + emp1 +
                  emp2,data=dat,family=binomial(link="probit"))
# wald and LR tests
require(aod)
names(coef(m1)) # the 6-7 is employment
wald.test(b = coef(m1), Sigma = vcov(m1), Terms = 6:7)

m2<- glm(safe ~ religious + minority  + female + age,data=dat,family=binomial)
lr.test(m2,m1)

require(car)
linearHypothesis(m1, c("emp1 = emp2"))

m3 <- glm(safe ~ religious + minority  + female + age + I(emp1 +emp2),data=dat,family=binomial)
lr.test(m1,m3)

# odds ratio
list.coef(m1)
# Figure 5.3
ar <-list.coef(m1)$out
ar <- ar[-1,]
ggplot(ar,aes(y=reorder(variables,exp.b),x=exp.b,xmax=ul.exp.b,xmin=ll.exp.b)) + theme_bw() +
  geom_pointrange() + labs(x="Odds Ratio (95% CI)",y="") +
  geom_vline(xintercept=1,linetype=2,color="grey50") +
  scale_y_discrete(labels=c("female"="Female",
                            "minority"="Minority",
                            "religious"="Religious",
                            "age"="Age",
                            "emp1"="Employed",
                            "emp2"="Self-Employed"))
p1<-ggplot(ar,aes(y=reorder(variables,exp.b),x=exp.b,xmax=ul.exp.b,xmin=ll.exp.b)) + theme_bw() +
  geom_pointrange() + labs(x="Odds Ratio (95% CI)",y="") +
  geom_vline(xintercept=1,linetype=2,color="grey50") +
  scale_y_discrete(labels=c("female"="Female",
                            "minority"="Minority",
                            "religious"="Religious",
                            "age"="Age",
                            "emp1"="Employed",
                            "emp2"="Self-Employed"))
#ggsave("Fig5.3.eps", plot = p1, device=cairo_ps)


# predicted probabilities - gender
design <- margins.des(m1,ivs=expand.grid(female=c(0,1)))
pdat <- margins.dat(m1,design)
pdat

first.diff.fitted(m1, design, compare=c(1,2))

# Using some dplyr (tidyverse) to make this better for plotting
pdat<-mutate(pdat,xaxs=c("Male","Female"))
pdat$xaxs<-factor(pdat$xaxs,levels=c("Male","Female"))


ggplot(pdat,aes(x=xaxs,y=fitted,ymin=ll,ymax=ul)) + theme_bw() +
  geom_pointrange() + labs(x="",y="Pr(Feels Safe)") +
  scale_y_continuous(limits=c(.6,.9))


# example 2, age
design <- margins.des(m1,expand_grid(age=seq(20,80,5)))
pdat <- margins.dat(m1,design)
# Figure 5.4
ggplot(pdat,aes(x=age,y=fitted,ymax=ul,ymin=ll)) + theme_bw() +
  geom_line() + geom_ribbon(alpha=.2) + labs(x="Age",y="Pr(Feels Safe)") +
  scale_y_continuous(limits=c(.7,.9))
#ggsave("Fig5.4.eps", plot = ggplot(pdat,aes(x=age,y=fitted,ymax=ul,ymin=ll)) + theme_bw() +
         geom_line() + geom_ribbon(alpha=.2) + labs(x="Age",y="Pr(Feels Safe)") +
         scale_y_continuous(limits=c(.7,.9)), device=cairo_ps)


first.diff.fitted(m1, design, compare=c(3,10))



require(marginaleffects)
summary(marginaleffects(m1))
# or, require(margins); summary(margins(m1)) # the latter is reported



### Diagnostics
diags <- diagn(m1)
diags[1:10,c(1,3,6)]

# residual plot; Figure 5.5
ggplot(diags,aes(x=obs,y=devres)) + geom_point() +
  geom_hline(yintercept=0,color="red") + xlab("Observation") +
  ylab("Deviance Residual") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme_bw()

dat.5 <- data.frame(dat,diags)
m1.sr <- glm(safe ~ religious + minority  + female + age + emp1 + emp2,data=filter(dat.5,devres >= -2 & devres <=2),family="binomial")
round(cbind(full.coef=coef(m1),sr.coef=coef(m1.sr),full.se=sqrt(diag(vcov(m1))),sr.se=sqrt(diag(vcov(m1.sr))),
            full.z=coef(m1)/sqrt(diag(vcov(m1))),sr.z=coef(m1.sr)/sqrt(diag(vcov(m1.sr)))),3)

round(diags[1:10,"deltabeta"],4)

# Figure 5.6
ggplot(diags,aes(x=obs,y=deltabeta)) + geom_point() +
  geom_hline(yintercept=4/nrow(diagn),color="red") +
  geom_text(data=filter(diags,deltabeta>(4/nrow(diags))),aes(x=obs,y=deltabeta,label=obs),vjust=-.6) +
  xlab("Observation") + ylab("Cook's D") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme_bw() +
  scale_y_continuous()




logLik(m1) # LL #1
m0<-glm(safe ~ 1,data=dat,family=binomial)
logLik(m0) # LL #2
deviance(m1)
lr.test(m1,m0)
require("pscl")
pR2(m1) # Fewer than Stata, but none are useful...
AIC(m1)
BIC(m1)

AIC(m2,m1)
BIC(m2,m1)

