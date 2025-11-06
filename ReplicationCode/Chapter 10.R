###
# Updated 6/17/2024 following release of catregs on CRAN
###

rm(list=ls())
require(tidyverse)
# install.packages("catregs")
require(catregs)
data("LF06travel")
LF06travel <- mutate(LF06travel,mode2=mode)
LF06travel$mode2[which(LF06travel$mode==1)]<-"Train"
LF06travel$mode2[which(LF06travel$mode==2)]<-"Bus"
LF06travel$mode2[which(LF06travel$mode==3)]<-"Car"

table(filter(LF06travel,choice==1)$mode2 )
pdat <- data.frame(y=table(filter(LF06travel,choice==1)$mode2 ))
ggplot(pdat,aes(x=y.Var1,y=y.Freq)) + theme_bw() + geom_col()

require(Epi)
m1 <- clogistic(choice ~ car + bus + time , strata=id, data=LF06travel)
m1

design <- data.frame(car=c(0,0,1),bus=c(0,1,0),time=c(60,60,60))
design
ma1 <- margins.dat.clogit(m1,design)
ma1
ma1<-mutate(ma1,type=c("Train","Bus","Car"))
# Figure 10.1
ggplot(ma1,aes(x=type,y=probs,ymin=ll,ymax=ul)) + theme_bw() +
  geom_point() + geom_errorbar(width=.1) + labs(x="",y="Probabliity of Selection")

compare.margins(margins=c(ma1$probs)[1:2],margins.ses=c(ma1$se)[1:2])
compare.margins(margins=c(ma1$probs)[c(1,3)],margins.ses=c(ma1$se)[c(1,3)])
compare.margins(margins=c(ma1$probs)[2:3],margins.ses=c(ma1$se)[2:3])







data("gss2016")
dim(gss2016)
sum(table(filter(gss2016,occ==1)$sclass))
round(table(filter(gss2016,occ==1)$sclass)/sum(table(filter(gss2016,occ==1)$sclass)),3)
m2 <- clogistic(occ ~ skmanual+ selfemp + service + proflow + profhigh, strata=id, data=gss2016)
m2

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
margins.dat.clogit(m2,design1)

m3 <- clogistic(occ ~ skmanual*pclass + selfemp*pclass + service*pclass + proflow*pclass + profhigh*pclass, strata=id, data=gss2016)
m3

# Reference category for father's occupation is farm
design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0)
margins.dat.clogit(m3,design1)
mar1<-margins.dat.clogit(m3,design1,rounded=3)
mar1


design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*1,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*1,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*1,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*1,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*1,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0)

mar2<-margins.dat.clogit(m3,design1,rounded=3)
mar2

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*1,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*1,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*1,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*1,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*1,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0)

mar3<-margins.dat.clogit(m3,design1,rounded=3)
mar3

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*1,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*1,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*1,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*1,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*1,up.lp=profhigh*0,up.up=profhigh*0)

mar4<-margins.dat.clogit(m3,design1,rounded=3)
mar4

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*1,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*1,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*1,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*1,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*1,up.up=profhigh*0)

mar5<-margins.dat.clogit(m3,design1,rounded=3)
mar5

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*1,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*1,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*1,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*1,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*1)

mar6<-margins.dat.clogit(m3,design1,rounded=3)
mar6



margins <- rbind(mar1,mar2,mar3,mar4,mar5,mar6)
margins <- mutate(margins,occs=rep(c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"),6),
                  foccs=rep(c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"),each=6))
margins$occs<-factor(margins$occs,levels=c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"))
margins$foccs<-factor(margins$foccs,levels=c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"))
# Figure 10.2
ggplot(margins,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw() + theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange") +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability") +
  geom_text(aes(label = probs))

# compare - 36,5
compare.margins(margins=margins$probs[c(36,6)],margins.ses=margins$se[c(36,6)])


m4 <- clogistic(occ ~ skmanual*pclass + selfemp*pclass + service*pclass + proflow*pclass + profhigh*pclass + skmanual*educ + selfemp*educ + service*educ + proflow*educ + profhigh*educ, strata=id, data=gss2016)
m4

# Margins when education is 12
design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*12,
                  se.e=selfemp*12,
                  sr.e=service*12,
                  lp.e=proflow*12,
                  up.e=profhigh*12)
margins.dat.clogit(m4,design1)
mar1<-margins.dat.clogit(m4,design1)
mar1

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*1,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*1,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*1,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*1,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*1,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*12,
                  se.e=selfemp*12,
                  sr.e=service*12,
                  lp.e=proflow*12,
                  up.e=profhigh*12)
mar2<-margins.dat.clogit(m4,design1)
mar2

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*1,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*1,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*1,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*1,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*1,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*12,
                  se.e=selfemp*12,
                  sr.e=service*12,
                  lp.e=proflow*12,
                  up.e=profhigh*12)
mar3<-margins.dat.clogit(m4,design1)
mar3

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*1,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*1,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*1,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*1,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*1,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*12,
                  se.e=selfemp*12,
                  sr.e=service*12,
                  lp.e=proflow*12,
                  up.e=profhigh*12)
mar4<-margins.dat.clogit(m4,design1)
mar4



design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*1,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*1,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*1,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*1,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*1,up.up=profhigh*0,
                  sk.e=skmanual*12,
                  se.e=selfemp*12,
                  sr.e=service*12,
                  lp.e=proflow*12,
                  up.e=profhigh*12)
mar5<-margins.dat.clogit(m4,design1)
mar5


design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*1,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*1,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*1,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*1,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*1,
                  sk.e=skmanual*12,
                  se.e=selfemp*12,
                  sr.e=service*12,
                  lp.e=proflow*12,
                  up.e=profhigh*12)
margins.dat.clogit(m4,design1)
mar6<-margins.dat.clogit(m4,design1)
mar6

margins <- rbind(mar1,mar2,mar3,mar4,mar5,mar6)
margins <- mutate(margins,occs=rep(c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"),6),
                  foccs=rep(c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"),each=6))
margins$occs<-factor(margins$occs,levels=c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"))
margins$foccs<-factor(margins$foccs,levels=c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"))
ggplot(margins,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw() + theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange",limits=c(0,.4)) +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability",title="Respondent has 12 years of Education") +
  geom_text(aes(label = probs))
p1.ed <- ggplot(margins,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw() +
  scale_fill_gradient(low="purple",high="orange",limits=c(0,.4)) +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability",title="Respondent has 12 years of Education") +
  geom_text(aes(label = probs))


# Margins when education is 16
design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*16,
                  se.e=selfemp*16,
                  sr.e=service*16,
                  lp.e=proflow*16,
                  up.e=profhigh*16)
mar1<-margins.dat.clogit(m4,design1)
mar1

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*1,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*1,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*1,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*1,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*1,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*16,
                  se.e=selfemp*16,
                  sr.e=service*16,
                  lp.e=proflow*16,
                  up.e=profhigh*16)
mar2<-margins.dat.clogit(m4,design1)
mar2

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*1,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*1,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*1,sr.sr=service*0,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*1,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*1,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*16,
                  se.e=selfemp*16,
                  sr.e=service*16,
                  lp.e=proflow*16,
                  up.e=profhigh*16)
mar3<-margins.dat.clogit(m4,design1)
mar3

design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*1,sk.lp=skmanual*0,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*1,se.lp=selfemp*0,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*1,sr.lp=service*0,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*1,lp.lp=proflow*0,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*1,up.lp=profhigh*0,up.up=profhigh*0,
                  sk.e=skmanual*16,
                  se.e=selfemp*16,
                  sr.e=service*16,
                  lp.e=proflow*16,
                  up.e=profhigh*16)
mar4<-margins.dat.clogit(m4,design1)
mar4



design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*1,sk.up=skmanual*0,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*1,se.up=selfemp*0,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*1,sr.up=service*0,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*1,lp.up=proflow*0,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*1,up.up=profhigh*0,
                  sk.e=skmanual*16,
                  se.e=selfemp*16,
                  sr.e=service*16,
                  lp.e=proflow*16,
                  up.e=profhigh*16)
mar5<-margins.dat.clogit(m4,design1)
mar5


design1 <- data.frame(skmanual=c(0,1,0,0,0,0),selfemp=c(0,0,1,0,0,0),service=c(0,0,0,1,0,0),proflow=c(0,0,0,0,1,0),profhigh=c(0,0,0,0,0,1))
design1 <- mutate(design1,
                  sk.sk=skmanual*0,sk.se=skmanual*0,sk.sr=skmanual*0,sk.lp=skmanual*0,sk.up=skmanual*1,
                  se.sk=selfemp*0,se.se=selfemp*0,se.sr=selfemp*0,se.lp=selfemp*0,se.up=selfemp*1,
                  sr.sk=service*0,sr.se=service*0,sr.sr=service*0,sr.lp=service*0,sr.up=service*1,
                  lp.sk=proflow*0,lp.se=proflow*0,lp.sr=proflow*0,lp.lp=proflow*0,lp.up=proflow*1,
                  up.sk=profhigh*0,up.se=profhigh*0,up.sr=profhigh*0,up.lp=profhigh*0,up.up=profhigh*1,
                  sk.e=skmanual*16,
                  se.e=selfemp*16,
                  sr.e=service*16,
                  lp.e=proflow*16,
                  up.e=profhigh*16)
mar6<-margins.dat.clogit(m4,design1)
mar6

margins2 <- rbind(mar1,mar2,mar3,mar4,mar5,mar6)
margins2 <- mutate(margins2,occs=rep(c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"),6),
                  foccs=rep(c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"),each=6))
margins2$occs<-factor(margins2$occs,levels=c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"))
margins2$foccs<-factor(margins2$foccs,levels=c("Unskilled Manual","Skilled Manual","Self-Employed","Non-Manual/Service","Professional, Lower","Professional, Higher"))
ggplot(margins2,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw() + theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange",limits=c(0,.4)) +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability",title="Respondent has a High School Education") +
  geom_text(aes(label = probs))
p2.ed <- ggplot(margins2,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw()  +
  scale_fill_gradient(low="purple",high="orange",limits=c(0,.4)) +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability",title="Respondent has 16 years of Education") +
  geom_text(aes(label = probs))
require(ggpubr)
ggarrange(p1.ed,p2.ed,nrow=2,labels=c("A","B")) # Figure 10.3

compare.margins(margins=c(margins2$probs[6],margins$probs[6]),
                margins.ses=c(margins2$se[6],margins$se[6]))

compare.margins(margins=c(margins2$probs[15],margins$probs[15]),
                margins.ses=c(margins2$se[15],margins$se[15]))



###
### Rank-Order Logit. Only with No Ties.
rm(list=ls())

require(mlogit)
data("Game")
G <- dfidx(Game, varying = 1:12, choice = "ch", ranked = TRUE, idnames = c("chid", "alt"))
m1<- mlogit(ch ~ own | hours + age, G, reflevel = "PC")
summary(m1)
require(catregs)

design <- data.frame(own=c(0,0,0,0,0,0),hours=0,age=35)
margins.dat.clogit(m1,design)
design2 <- data.frame(own=c(0,0,0,1,0,0),hours=0,age=35)
margins.dat.clogit(m1,design2)
design3 <- data.frame(own=c(1,0,0,1,0,1),hours=0,age=35)
margins.dat.clogit(m1,design3)


mean(Game$own.Xbox)
mean(Game$own.PlayStation)
mean(Game$own.PSPortable)
mean(Game$own.GameCube)
mean(Game$own.GameBoy)
mean(Game$own.PC)
mean(Game$hours)
mean(Game$age)
sd(Game$hours)

design <- data.frame(own=c(.88,.13,.08,.32,.1,.13),hours=3.88,age=22.23)
mar1<-margins.dat.clogit(m1,design)
design2 <- data.frame(own=c(.88,.13,.08,.32,.1,.13),hours=3.88 + 5.01,age=22.23)
mar2<-margins.dat.clogit(m1,design2)
mar1<- mutate(mar1,type=c("Computer","GameBoy","GameCube","PlayStation","PSPortable","Xbox"))
mar2<- mutate(mar2,type=c("Computer","GameBoy","GameCube","PlayStation","PSPortable","Xbox"))
p1<-ggplot(mar1,aes(x=type,y=probs,ymin=ll,ymax=ul)) + theme_bw() +
  geom_point() + geom_errorbar(width=.1) + labs(x="",y="Probability of Selection",title="Average Hours Playing") +
  scale_y_continuous(limits=c(0,.65))
p2<-ggplot(mar2,aes(x=type,y=probs,ymin=ll,ymax=ul)) + theme_bw() +
  geom_point() + geom_errorbar(width=.1) + labs(x="",y="Probability of Selection",title="Average + Std. Dev. Hours Playing") +
  scale_y_continuous(limits=c(0,.65))
require(ggpubr)
p4<-ggarrange(p1,p2,labels=c("A","B"))
p4 # Figure 10.4

# Computer
compare.margins(margins=c(mar2$probs[1],mar1$probs[1]),
                margins.ses=c(mar2$se[1],mar1$se[1]))
# Gameboy
compare.margins(margins=c(mar2$probs[2],mar1$probs[2]),
                margins.ses=c(mar2$se[2],mar1$se[2]))
# Gamecube
compare.margins(margins=c(mar2$probs[3],mar1$probs[3]),
                margins.ses=c(mar2$se[3],mar1$se[3]))
# Playstation
compare.margins(margins=c(mar2$probs[4],mar1$probs[4]),
                margins.ses=c(mar2$se[4],mar1$se[4]))
# PSPortable
compare.margins(margins=c(mar2$probs[5],mar1$probs[5]),
                margins.ses=c(mar2$se[5],mar1$se[5]))
# Xbox
compare.margins(margins=c(mar2$probs[6],mar1$probs[6]),
                margins.ses=c(mar2$se[6],mar1$se[6]))


###
# Bootstrap uncertainty
###
m1<- mlogit(ch ~ own | hours + age, G, reflevel = "PC")
summary(m1)
design <- data.frame(own=c(0,0,0,0,0,0),hours=0,age=35)
obs.predictions<-as.matrix(predict(m1,newdata=design),nc=1)


bs.dist <- matrix(NA,nr=1000,nc=length(obs.predictions))
for(i in 1:1000){
  set.seed(i); Game.bs <- Game[sample(1:nrow(Game),.9*nrow(Game),replace=TRUE),] # code in the bracket sample 90% of respondents with replacement
  G.bs <- dfidx(Game.bs, varying = 1:12, choice = "ch", ranked = TRUE, idnames = c("chid", "alt"))
  m1.bs<- mlogit(ch ~ own | hours + age, G.bs, reflevel = "PC")
  bs.dist
  bs.dist[i,]<-predict(m1.bs,newdata=design)}
bs.dist<-apply(bs.dist,2,FUN="sort")

obs.predictions
pdat <- data.frame(observed=obs.predictions,
                   vars=rownames(obs.predictions),
                   mins=bs.dist[.025*nrow(bs.dist),],
                   maxs=bs.dist[.975*nrow(bs.dist),])
ggplot(pdat,aes(x=vars,y=observed,ymin=mins,ymax=maxs)) +
  theme_bw() + geom_pointrange() + labs(x="",y="Probability of Selection")










###
### Another conditional logistic regression example. The data are too old for the book
###

# Conditional Logistic Regression
rm(list=ls())
require(catregs)
data(logan)
head(logan2)
require(Epi)
require(tidyverse)
filter(logan2,id==1)

table(filter(logan2,case==1)$occupation)
table(filter(logan2,case==1)$occupation)/sum(table(filter(logan2,case==1)$occupation))

logan2 <- mutate(logan2,sales=ifelse(operatives==0 & craftsmen==0 & farm==0 & professional==0,1,0))
m5 <- clogistic(case ~ operatives + craftsmen + sales + professional  , strata=id, data=logan2)
design1 <- data.frame(operatives=c(0,1,0,0,0),craftsmen=c(0,0,1,0,0),sales=c(0,0,0,1,0),professional=c(0,0,0,0,1))
margins.dat.clogit(m5,design1)
# Baseline probability of being in any occupational category; same as the proportions in the table statement...

m6 <- clogistic(case ~ operatives*focc + craftsmen*focc + sales*focc + professional*focc  , strata=id, data=logan2)
m6

# Reference category for father's occupation is farm
design1 <- data.frame(operatives=c(0,1,0,0,0),craftsmen=c(0,0,1,0,0),sales=c(0,0,0,1,0),professional=c(0,0,0,0,1))
design1 <- mutate(design1,oo=operatives*0,oc=operatives*0,os=operatives*0,op=operatives*0,
                  co=craftsmen*0,cc=craftsmen*0,cs=craftsmen*0,cp=craftsmen*0,
                  so=sales*0,sc=sales*0,ss=sales*0,sp=sales*0,
                  po=professional*0,pc=professional*0,ps=professional*0,pp=professional*0)
mar1<-margins.dat.clogit(m6,design1,rounded=2)
mar1

design1 <- data.frame(operatives=c(0,1,0,0,0),craftsmen=c(0,0,1,0,0),sales=c(0,0,0,1,0),professional=c(0,0,0,0,1))
design1 <- mutate(design1,oo=operatives*1,oc=operatives*0,os=operatives*0,op=operatives*0,
                  co=craftsmen*1,cc=craftsmen*0,cs=craftsmen*0,cp=craftsmen*0,
                  so=sales*1,sc=sales*0,ss=sales*0,sp=sales*0,
                  po=professional*1,pc=professional*0,ps=professional*0,pp=professional*0)
mar2<-margins.dat.clogit(m6,design1,rounded=2)
mar2

design1 <- data.frame(operatives=c(0,1,0,0,0),craftsmen=c(0,0,1,0,0),sales=c(0,0,0,1,0),professional=c(0,0,0,0,1))
design1 <- mutate(design1,oo=operatives*0,oc=operatives*1,os=operatives*0,op=operatives*0,
                  co=craftsmen*0,cc=craftsmen*1,cs=craftsmen*0,cp=craftsmen*0,
                  so=sales*0,sc=sales*1,ss=sales*0,sp=sales*0,
                  po=professional*0,pc=professional*1,ps=professional*0,pp=professional*0)
mar3<-margins.dat.clogit(m6,design1,rounded=2)
mar3

design1 <- data.frame(operatives=c(0,1,0,0,0),craftsmen=c(0,0,1,0,0),sales=c(0,0,0,1,0),professional=c(0,0,0,0,1))
design1 <- mutate(design1,oo=operatives*0,oc=operatives*0,os=operatives*1,op=operatives*0,
                  co=craftsmen*0,cc=craftsmen*0,cs=craftsmen*1,cp=craftsmen*0,
                  so=sales*0,sc=sales*0,ss=sales*1,sp=sales*0,
                  po=professional*0,pc=professional*0,ps=professional*1,pp=professional*0)
mar4<-margins.dat.clogit(m6,design1,rounded=2)
mar4

design1 <- data.frame(operatives=c(0,1,0,0,0),craftsmen=c(0,0,1,0,0),sales=c(0,0,0,1,0),professional=c(0,0,0,0,1))
design1 <- mutate(design1,oo=operatives*0,oc=operatives*0,os=operatives*0,op=operatives*1,
                  co=craftsmen*0,cc=craftsmen*0,cs=craftsmen*0,cp=craftsmen*1,
                  so=sales*0,sc=sales*0,ss=sales*0,sp=sales*1,
                  po=professional*0,pc=professional*0,ps=professional*0,pp=professional*1)
mar5<-margins.dat.clogit(m6,design1,rounded=2)
mar5

# mar1 is when dad is farmer
# mar2 is when dad is operative
# mar3 is when dad is craftsmen
# mar4 is when dad is sales
# mar5 is when dad is professional
margins <- rbind(mar1,mar2,mar3,mar4,mar5)
margins <- mutate(margins,occs=rep(c("Farmer","Operatives","Craft","Sales","Professional"),5),
                  foccs=rep(c("Farmer","Operatives","Craft","Sales","Professional"),each=5))
ggplot(margins,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw() + theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange") +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability") +
  geom_text(aes(label = probs))
p4<-ggplot(margins,aes(x=occs,y=foccs,fill=probs)) + geom_tile() +
  theme_bw() + theme(legend.position="bottom") +
  scale_fill_gradient(low="purple",high="orange") +
  labs(x="Respondent's Occupational Category",y="Father's Occupational Category",
       fill="Probability") +
  geom_text(aes(label = probs))
#ggsave("Fig10.4.eps", plot = p4)

# Compare probability of being professional by prof/farm father
compare.margins(margins=margins$probs[c(25,5)],margins.ses=margins$se[c(25,5)])


