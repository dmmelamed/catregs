###
# Updated 6/17/2024 following release of catregs on CRAN
###


rm(list=ls())
require(tidyverse)
# install.packages("catregs)
require(catregs)
data(essUK)

X <- filter(essUK,country=="United Kingdom")
table(X$walk.alone.dark)

X <- mutate(X,safe = ifelse(walk.alone.dark=="Safe" | walk.alone.dark==
                              "Very safe",1,0),cated=NA)
X$cated[which(X$education<13)]<-"HS Or Less"
X$cated[which(X$education > 12 & X$education < 16)]<-"Some College"
X$cated[which(X$education == 16)]<-"BA/BS"
X$cated[which(X$education > 16)]<-"Grad School"
X$cated <- factor(X$cated,levels=c("HS Or Less","Some College","BA/BS","Grad School"))


###
# a.	Univariate
# i.	Descriptive Statistics
###
# Table 4.1
table(X$gender)
table(X$gender)/nrow(X)


table(X$cated)
table(X$cated)/sum(table(X$cated))

X %>% drop_na(cated) %>% group_by(cated) %>%
  summarize(n=n()) %>%
  mutate(prop=n/sum(n),
         cum.prop=cumsum(prop))

###
# # b.	Bivariate
# i.	Cross-tabulation
# ii.	Relative Risk and the Odds Ratio
# iii.	Bivariate Statistical Tests
###


# Table 4.2
table(X$gender,X$safe) # Cross-tabulations

# Table 4.3
table(X$gender,X$safe)/sum(table(X$gender,X$safe)) # Proportions
rowSums(table(X$gender,X$safe)/sum(table(X$gender,X$safe)))
colSums(table(X$gender,X$safe)/sum(table(X$gender,X$safe)))
table(X$gender,X$safe)/rowSums(table(X$gender,X$safe)) # Conditional Proportions (Conditional on the row variable)
# Difference in proportions
props<- (table(X$gender,X$safe)/rowSums(table(X$gender,X$safe)))[,2] # Conditional Proportions (Conditional on the row variable)
props[2] - props[1]


# Relative risk, ratio of proportions
props[1]/props[2]
props[2]/props[1]


# Odds
props[1]/(1-props[1]) # Women are 2.03 times more likely to report feeling safe than not safe
props[2]/(1-props[2]) # Men are 5.95 times more likely to report feeling safe than not safe


# Odds Ratio
o1<-props[1]/(1-props[1]) # For women
o2<-props[2]/(1-props[2]) # For men
o1/o2 # Women are .34 times as likely as men to report feeling safe walking alone at night.
# Women are 66% less likely to report...
o2/o1 # Men are 2.93 times more likely than women to report feeling safe walking alone at night
1/(o2/o1) # Same as o1/o2





# iii.	Bivariate Statistical Tests
# Chi-squared here. Fit versus not fit.

# Start with 2 x 2 table
table(X$gender,X$safe)
chisq.test(table(X$gender,X$safe))
chisq.test(table(X$gender,X$safe),correct=FALSE)


c1<-chisq.test(table(X$gender,X$safe))
c1
names(c1)
c1$residuals
c1$observed
c1$expected
(abs(c1$observed - c1$expected)-.5)^2/c1$expected # Chi-Squared Components
# Interpreting this is weird since it's a 2x2. Will do so below.
sum((abs(c1$observed - c1$expected)-.5)^2/c1$expected)
c2<-chisq.test(table(X$gender,X$safe),correct=FALSE)
(c2$observed - c2$expected)^2/c2$expected # Chi-Squared Components
sum((c2$observed - c2$expected)^2/c2$expected)
pchisq(100.87,1,lower.tail=FALSE)



# Bigger table



table(X$can.trust.people)
X <-mutate(X,trust2=NA)
X$trust2[which(X$can.trust.people<5)] <-1
X$trust2[which(X$can.trust.people>5 & X$can.trust.people<7)] <-2
X$trust2[which(X$can.trust.people>6)] <-3
X <-mutate(X,life.sat2=NA)
X$life.sat2[which(X$life.satisfaction<5)] <-1
X$life.sat2[which(X$life.satisfaction==5 | X$life.satisfaction==6)] <-2
X$life.sat2[which(X$life.satisfaction==7)] <-3
X$life.sat2[which(X$life.satisfaction==8)] <-4
X$life.sat2[which(X$life.satisfaction==9 | X$life.satisfaction==10)] <-5
table(X$life.sat2,X$trust2)
chisq.test(table(X$life.sat2,X$trust2))

c3 <- chisq.test(table(X$life.sat2,X$trust2))

res1<-(c3$observed - c3$expected)^2/c3$expected # Chi-Squared Components
res1 <- data.frame(res1)

which(data.frame(c3$residuals)$Freq < 0)
res1$Freq[which(data.frame(c3$residuals)$Freq < 0)] <- res1$Freq[which(data.frame(c3$residuals)$Freq < 0)] * -1
# Figure 4.1
ggplot(res1,aes(x=Var2,y=Var1,fill=Freq,label=round(Freq,1))) + geom_tile() + theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange") +
  labs(x="Trust in Others",y="Satisfied with Life",fill="Residual") +
  geom_text()  +
  scale_x_discrete(breaks=1:3,labels=c("Few People","Some People","Most People")) +
  scale_y_discrete(breaks=1:5,labels=c("Not Satisfied",
                                         "A Little Satisfied",
                                         "Somewhat Satisfied",
                                         "Very Satisfied",
                                         "Completely Satisfied"))
p1<- ggplot(res1,aes(x=Var2,y=Var1,fill=Freq,label=round(Freq,1))) + geom_tile() + theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange") +
  labs(x="Trust in Others",y="Satisfied with Life",fill="Residual") +
  geom_text()  +
  scale_x_discrete(breaks=1:3,labels=c("Few People","Some People","Most People")) +
  scale_y_discrete(breaks=1:5,labels=c("Not Satisfied",
                                       "A Little Satisfied",
                                       "Somewhat Satisfied",
                                       "Very Satisfied",
                                       "Completely Satisfied"))
#ggsave("Fig4.1.eps", plot = p1)


# c.	Multivariate: Log-linear models
# i.	The Log-Link
# ii.	Model selection and interpretation
# iii.	Model Comparison with MLE

# Start w loglin for 2-way table. Table 4.6
table(X$gender,X$safe)
lldat <- data.frame(table(X$gender,X$safe))
colnames(lldat)<-c("Sex","Safe","Count")
ll.m1 <- glm(Count ~ Sex*Safe,data=lldat,family="poisson")
summary(ll.m1)
round(predict(ll.m1,type="response"),1)
ll.m2 <- glm(Count ~ Sex + Safe,data=lldat,family="poisson")
summary(ll.m2)
round(predict(ll.m1,type="response"),1)
require(stargazer) # Load the package
stargazer(ll.m1,ll.m2, type = "html", out="loglins.html",
          star.cutoffs=c(.05,.01,.001),omit.stat=c("aic","n","bic"))


AIC(ll.m2)
-2*logLik(ll.m2) + 6

BIC(ll.m2)
-2*logLik(ll.m2) + log(nrow(ll.m2$model))*3
anova(ll.m1,ll.m2)
#Compare to traditional Chi-squared statistic
obs <-table(X$gender,X$safe)
exp <- outer(rowSums(obs),colSums(obs))/sum(obs)
sum(((obs-exp)^2)/exp)
pchisq(sum(((obs-exp)^2)/exp),1,lower.tail=FALSE)

# 3-way table. Table 4.7
table(X$cated,X$gender,X$safe)
lldat2 <- data.frame(table(X$cated,X$gender,X$safe))
colnames(lldat2)[1:3] <- c("E","F","S")

# Table 4.8
ll.sat <- glm(Freq ~ E*F*S,data=lldat2,family="poisson")
ll.main <- glm(Freq ~ E*S + F*S,data=lldat2,family="poisson")
anova(ll.sat,ll.main)
ll.no.ed <- glm(Freq ~ E + F*S,data=lldat2,family="poisson")
anova(ll.main,ll.no.ed)
ll.no.sex <- glm(Freq ~ E*S + F,data=lldat2,family="poisson")
anova(ll.main,ll.no.sex)

# Marginal Proportions, Table 4.9
predict(ll.main,type="response")
lldat2 <- mutate(lldat2,yhat=predict(ll.main,type="response"))
lldat2[1,5] / (lldat2[1,5] + lldat2[9,5] )
lldat2[9,5] / (lldat2[1,5] + lldat2[9,5] )
lldat2[2,5] / (lldat2[2,5] + lldat2[10,5] )
lldat2[10,5] / (lldat2[10,5] + lldat2[2,5] )
lldat2[3,5] / (lldat2[3,5] + lldat2[11,5] )
lldat2[11,5] / (lldat2[11,5] + lldat2[3,5] )
lldat2[4,5] / (lldat2[4,5] + lldat2[12,5] )
lldat2[12,5] / (lldat2[12,5] + lldat2[4,5] )
lldat2[5,5] / (lldat2[5,5] + lldat2[13,5] )
lldat2[13,5] / (lldat2[13,5] + lldat2[5,5] )
lldat2[6,5] / (lldat2[6,5] + lldat2[14,5] )
lldat2[14,5] / (lldat2[14,5] + lldat2[6,5] )
lldat2[7,5] / (lldat2[7,5] + lldat2[15,5] )
lldat2[15,5] / (lldat2[15,5] + lldat2[7,5] )
lldat2[8,5] / (lldat2[8,5] + lldat2[16,5] )
lldat2[16,5] / (lldat2[16,5] + lldat2[8,5] )

(.78/(1-.78)) / (.56/(1-.56))
(.91/(1-.91)) / (.79/(1-.79))


# 4-way table.. Table 4.10
table(X$cated,X$minority,X$gender,X$safe)
lldat3 <- data.frame(table(X$cated,X$minority,X$gender,X$safe))
colnames(lldat3)[1:4] <- c("E","M","F","S")
# Saturated
ll.sat <- glm(Freq ~ E*M*F*S,data=lldat3,family="poisson")
length(coef(ll.sat))
logLik(ll.sat)
BIC(ll.sat);AIC(ll.sat)

ll.2 <- glm(Freq ~ M*F*S  + E*F*S + E*M*S + E*M*F,data=lldat3,family="poisson")
anova(ll.sat,ll.2) # Can constrain 4-way interaction
length(coef(ll.2))
logLik(ll.2)
BIC(ll.2);AIC(ll.2)

ll.3 <- glm(Freq ~ E*F*S + E*M*S + E*M*F,data=lldat3,family="poisson")
length(coef(ll.3))
logLik(ll.3)
BIC(ll.3);AIC(ll.3)

ll.4 <- glm(Freq ~ M*F*S   + E*M*S + E*M*F,data=lldat3,family="poisson")
length(coef(ll.4))
logLik(ll.4)
BIC(ll.4);AIC(ll.4)

ll.5 <- glm(Freq ~ M*F*S  + E*F*S + E*M*F,data=lldat3,family="poisson")
length(coef(ll.5))
logLik(ll.5)
BIC(ll.5);AIC(ll.5)

ll.6 <- glm(Freq ~ M*F*S  + E*F*S + E*M*S,data=lldat3,family="poisson")
length(coef(ll.6))
logLik(ll.6)
BIC(ll.6);AIC(ll.6)

anova(ll.2,ll.3)
anova(ll.2,ll.4)
anova(ll.2,ll.5)
anova(ll.2,ll.6)

ll.7 <- glm(Freq ~ M*F*S  + E*F + E*M + E*S,data=lldat3,family="poisson")
length(coef(ll.7))
logLik(ll.7)
BIC(ll.7);AIC(ll.7)
anova(ll.2,ll.7)
ll.8 <- glm(Freq ~ M*F*S  +  E*M + E*S,data=lldat3,family="poisson")
length(coef(ll.8))
logLik(ll.8)
BIC(ll.8);AIC(ll.8)
anova(ll.7,ll.8) # Sig
ll.9 <- glm(Freq ~ M*F*S  + E*F +  E*S,data=lldat3,family="poisson")
length(coef(ll.9))
logLik(ll.9)
BIC(ll.9);AIC(ll.9)
anova(ll.7,ll.9)
ll.10 <- glm(Freq ~ M*F*S  + E*F + E*M ,data=lldat3,family="poisson")
length(coef(ll.10))
logLik(ll.10)
BIC(ll.10);AIC(ll.10)
anova(ll.7,ll.10)

lldat3 <- mutate(lldat3,yhat=round(predict(ll.7,type="response"),2))
colnames(lldat3) <- c("Education","Minority","Female","Safe","N","Model 7 Fitted Values")
lldat3

write.csv(lldat3,"ch4dat.csv") # Appendix A


