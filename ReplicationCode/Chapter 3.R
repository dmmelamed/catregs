
rm(list=ls())
require(tidyverse)
require(catregs)
data(essUK)

head(essUK)


ggplot(essUK,aes(x=can.trust.people)) + theme_bw() +
  geom_bar() + scale_x_continuous(breaks=seq(0,10,2)) +
  labs(x="Generalized Trust")

ggplot(essUK,aes(x=education)) + theme_bw() +
  geom_bar() + scale_x_continuous(breaks=seq(0,20,2),limits=c(0,20)) +
  labs(x="Education")

ggplot(essUK,aes(x=education,y=can.trust.people)) + theme_bw() +
  geom_jitter(alpha=.7) + scale_x_continuous(breaks=seq(2,10,2),limits=c(2,20)) +
  labs(y="Generalized Trust",x="Education") + scale_y_continuous(breaks=seq(2,10,2)) +
  geom_smooth(method="lm",se=FALSE)

p1 <-ggplot(essUK,aes(x=can.trust.people)) + theme_bw() +
  geom_bar() + scale_x_continuous(breaks=seq(0,10,2)) +
  labs(x="Generalized Trust")

p2 <- ggplot(essUK,aes(x=education)) + theme_bw() +
  geom_bar() + scale_x_continuous(breaks=seq(0,20,2),limits=c(0,20)) +
  labs(x="Education")
require(ggpubr)
p1 <- ggarrange(p1,p2,labels=c("A","B"))
p3<- ggplot(essUK,aes(x=education,y=can.trust.people)) + theme_bw() +
  geom_jitter() + scale_x_continuous(breaks=seq(2,10,2),limits=c(2,20)) +
  labs(y="Generalized Trust",x="Education") + scale_y_continuous(breaks=seq(2,10,2)) +
  geom_smooth(method="lm",se=FALSE)
ggarrange(p1,p3,nrow=2,labels=c("","C")) # Figure 3.1
#p1<-ggarrange(p1,p3,nrow=2,labels=c("","C"))
#ggsave("Fig3.1.eps", plot = p1)

X <- select(essUK,can.trust.people,religious,minority,female,age,married,education)
X <- na.omit(X)

m1 <- lm(can.trust.people ~ education ,data=X)
m2 <- lm(can.trust.people ~ education + religious + age + minority + female + age,data=X)
summary(m1)
names(m1)
coef(m2)
sqrt(diag(vcov(m2)))
require(stargazer) # Load the package
stargazer(m1,m2, type = "html", out="ols1.html",
          star.cutoffs=c(.05,.01,.001),omit.stat=c("aic","n","bic")) # Table 3.1

require(broom)
out_m2 <- tidy(m2, conf.int=TRUE)
out_m2
ggplot(out_m2[-1,],aes(x=estimate,xmin=conf.low,xmax=conf.high,y=reorder(term,estimate))) + theme_bw() +
  geom_pointrange() + geom_vline(xintercept=0) +
  scale_y_discrete(labels=c("education"="Education","religious"="Religious","age"="Age","female"="Female","minority"="Minority")) +
  labs(y="",x="OLS Regression Coefficient (95% CI)") # Figure 3.2
p1 <-ggplot(out_m2[-1,],aes(x=estimate,xmin=conf.low,xmax=conf.high,y=reorder(term,estimate))) + theme_bw() +
  geom_pointrange() + geom_vline(xintercept=0) +
  scale_y_discrete(labels=c("education"="Education","religious"="Religious","age"="Age","female"="Female","minority"="Minority")) +
  labs(y="",x="OLS Regression Coefficient (95% CI)")
#ggsave("Fig3.2.eps", plot = p1)



design <- margins.des(m2,expand.grid(education=8:20))
round(design,2)
pdat <- margins.dat(m2,design)
ggplot(pdat,aes(x=education,y=fitted,ymin=ll,ymax=ul)) + theme_bw() +
  geom_line() + geom_ribbon(alpha=.2) + labs(x="Education",y="Predicted Generalized Trust") #Figure 3.3
p1<-ggplot(pdat,aes(x=education,y=fitted,ymin=ll,ymax=ul)) + theme_bw() +
  geom_line() + geom_ribbon(alpha=.2) + labs(x="Education",y="Predicted Generalized Trust")
#ggsave("Fig3.3.eps", plot = p1, device=cairo_ps)


# Table 3.2
m3 <- lm(can.trust.people ~ education*religious + age + minority + female + age,data=X)
m4 <- lm(can.trust.people ~ education*minority  + religious + age + minority + female + age,data=X)
require(stargazer) # Load the package
stargazer(m3,m4, type = "html", out="ols2.html",
          star.cutoffs=c(.05,.01,.001),omit.stat=c("aic","n","bic"))

table(X$religious)
mean(X$religious)
sd(X$religious)

# Figure 3.4
design2 <- margins.des(m3,expand.grid(education=8:20,religious=c(.54,3.6,6.66)))
pdat2 <- margins.dat(m3,design2)
pdat2$religious <- recode(pdat2$religious,".54"="Mean -1 SD","3.6"="Mean","6.66"="Mean +1 SD")
pdat2$religious <- factor(pdat2$religious,levels=c("Mean -1 SD","Mean","Mean +1 SD"))
ggplot(pdat2,aes(x=education,y=fitted,ymin=ll,ymax=ul,group=religious,fill=religious)) + theme_bw() + 
  geom_ribbon(alpha=.2) + geom_line() +
  theme(legend.position="bottom") +
  labs(x="Education",y="Predicted Generalized Trust",fill="Religiousness:") + scale_x_continuous(breaks=(seq(8,20,2)))
p1<-ggplot(pdat2,aes(x=education,y=fitted,ymin=ll,ymax=ul,group=religious,fill=religious)) + theme_bw() + 
  geom_ribbon(alpha=.2) + geom_line() +
  theme(legend.position="bottom") +
  labs(x="Education",y="Predicted Generalized Trust",fill="Religiousness:") + scale_x_continuous(breaks=(seq(8,20,2)))
#ggsave("Fig3.4.eps", plot = p1, device=cairo_ps)

design2.eg <- margins.des(m3,expand.grid(education=12:13,religious=3:4))
first.diff.fitted(m3,design2.eg,compare=c(2,1))
first.diff.fitted(m3,design2.eg,compare=c(4,3))
second.diff.fitted(m3,design2.eg,compare=c(4,3,2,1)) # The interaction effect is recovered







# Not in the book. Ilustrate model 2 in table 3.2
design3 <- margins.des(m4,expand.grid(education=2:20,minority=c(0,1)))
pdat3 <- margins.dat(m4,design3)
pdat3$minority <- recode(pdat3$minority,"0"="Not Minoritized","1"="Minoritized")
ggplot(pdat3,aes(x=education,y=fitted,ymin=ll,ymax=ul,group=minority,fill=minority)) + theme_bw() + 
  geom_ribbon(alpha=.2) + geom_line() +
  theme(legend.position="bottom") +
  labs(x="Education",y="Predicted Generalized Trust",fill="Minority Status:")



# Diagnostics
# our diagn function does not work in this context as there as good alternatives already
modfit <- augment(m2)
head(modfit) # The data, with the fitted.values and residuals appended
p1 <-ggplot(modfit,aes(x=.fitted, y = .std.resid)) + 
  geom_hline(yintercept=0,color="red") + geom_point()  + theme_bw() + labs(x="Predicted Value",y="Standardized Residual")
p1



pdat <- modfit %>% group_by(education) %>% summarize(means=mean(.std.resid),sds=sd(.std.resid))

p2<-ggplot(filter(pdat,education > 8 & education < 21),aes(x=education,y=means,ymin=means+sds,ymax=means-sds)) + theme_bw() +geom_hline(yintercept=0,color="red") +
  geom_pointrange() + labs(x="Education",y="Mean +/- SD of the Standardized Residual")
require(ggpubr)
p1<-ggarrange(p1,p2,labels=c("A","B"))
p1 # Figure 3.5
# ggsave("Fig3.5.eps", plot = p1, device=cairo_ps)



# Figure 3.6
contin_syn1 <- function(nobs = 200, xv =c (0,.6)){
  
  p<-length(xv)
  X<-cbind(matrix(rnorm(nobs * (p-1)),ncol = p-1))
  X56<-cbind(1,X)
  xb<-X56 %*% xv
  py<-rnorm(nobs,xb)
  out<-data.frame(cbind(py,X56[,-1]))
  names(out)<-c("py","x1")
  return(out)
  
}

dat <- contin_syn1(10000,xv=c(0,.5))
m1 <- lm(py~x1,data=dat)
sim.dat <- data.frame(dat,y=predict(m1,dat))
ggplot(sim.dat,aes(x=x1,y=y))+ theme_bw() +
  geom_smooth(method="lm",se=FALSE)+ 
  scale_y_continuous(limits=c(-.25,1.25),breaks=seq(-.2,1.2,.2)) +
  geom_hline(yintercept=0) + geom_hline(yintercept=1) +
  labs(y="Binary Response",x="Continuous Predictor")
p1 <- ggplot(sim.dat,aes(x=x1,y=y))+ theme_bw() +
  geom_smooth(method="lm",se=FALSE,color="gray55")+ 
  scale_y_continuous(limits=c(-.2,1.2),breaks=seq(-.2,1.2,.2)) +
  geom_hline(yintercept=0) + geom_hline(yintercept=1) +
  labs(y="Binary Response",x="Continuous Predictor") + scale_x_continuous(limits=c(-.1,2.1),breaks=c(0,.5,1,1.5,2),labels=c("-4","-2","0","2","4"))

dat <- mutate(dat,
              y2=as.numeric(py>0),
              p=glm(y2~x1,data=sim.dat,family="binomial")$fitted.values)
ggplot(dat,aes(x=x1,y=p))+ theme_bw() +
  geom_smooth(se=FALSE)+ 
  scale_y_continuous(limits=c(-.25,1.25),breaks=seq(-.2,1.2,.2)) +
  geom_hline(yintercept=0) + geom_hline(yintercept=1) +
  labs(y="Logit of the Binary Response",x="Continuous Predictor")
p2<-ggplot(dat,aes(x=x1,y=p))+ theme_bw() +
  geom_smooth(se=FALSE,color="gray55")+ 
  scale_y_continuous(limits=c(-.2,1.2),breaks=seq(-.2,1.2,.2)) +
  geom_hline(yintercept=0) + geom_hline(yintercept=1) + scale_x_continuous(limits=c(-4.2,4.2),breaks=seq(-4,4,2)) +
  labs(y="Logit of the Binary Response",x="Continuous Predictor")
require(ggpubr)
ggarrange(p1,p2,labels=c("A","B"))

p1 <- ggarrange(p1,p2,labels=c("A","B"))
#ggsave("Fig3.6.eps", plot = p1, device=cairo_ps)








