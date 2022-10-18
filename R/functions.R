lr.test<-function(full.model,reduced.model){
  # "Note: It does not matter if you correctly identified the Full and Reduced Models. The function detects it based on DF.")
  n.f <- length(predict(full.model))
  n.r <- length(predict(reduced.model))
  if(n.f != n.r){print("The models were not fit to the same data! WTF")}else{
    ll.f <-as.numeric(logLik(full.model))
    ll.r <- as.numeric( logLik(reduced.model))
    ll <- 2*abs(ll.r-ll.f)
    df <- abs(length(coef(full.model))-length(coef(reduced.model)))
    p.value <- round(pchisq(ll,df,lower.tail=FALSE),5)
    out <- data.frame("LL Full"=ll.f,"LL Reduced"=ll.r,
                      "G2/LR Statistic"=ll,"DF"=df,"p-value"=p.value)
    return(out)}}


list.coef<-function(model,rounded=3,alpha=.05){
  out<-matrix(0,nr=length(coef(model)),nc=10)
  out[,1]<-coef(model)
  vcov1<-vcov(model)[1:length(coef(model)),1:length(coef(model))]
  out[,2]<-sqrt(diag(vcov1))
  out[,3]<-out[,1]/out[,2]
  out[,4] <- out[,1]-qnorm(1-(alpha/2),lower.tail=TRUE)*out[,2]
  out[,5] <- out[,1]+qnorm(1-(alpha/2),lower.tail=TRUE)*out[,2]
  out[,6]<-dnorm(out[,3])
  out[,7]<-exp(coef(model))
  out[,8]<-exp(out[,4])
  out[,9]<-exp(out[,5])
  out[,10]<-100*(exp(coef(model))-1)
  out<-round(out,rounded)
  colnames(out)<-c("b","S.E.","Z","LL CI","UL CI","p-val","exp(b)","LL CI for exp(b)","UL CI for exp(b)","Percent")
  out <- data.frame(out,CI=paste(100*(1-alpha),"%"))
  rownames(out)<-names(coef(model))
  outp<-list(out=out)
  return(outp)}

compare.AMEs <- function(margins.matrix,seed=1234,rounded=3,nsim=10000){
  difference <- margins.matrix$AME[1] - margins.matrix$AME[2]
  if(difference>0){
    set.seed(seed); p.value<-sum((rnorm(nsim,mean=margins.matrix$AME[1],sd=margins.matrix$SE[1]) - rnorm(nsim,mean=margins.matrix$AME[2],sd=margins.matrix$SE[2])) < 0) /nsim
  }else{
    set.seed(seed); p.value<-sum((rnorm(nsim,mean=margins.matrix$AME[1],sd=margins.matrix$SE[1]) - rnorm(nsim,mean=margins.matrix$AME[2],sd=margins.matrix$SE[2])) > 0) /nsim
  }
  out <- list(difference=round(difference,rounded),p.value=round(p.value,rounded))
  return(out)
}


margins.des <- function(mod,ivs,excl="nonE"){
  X.mod <- mod$model[,-1]
  var.names<-colnames(X.mod)
  if(excl=="nonE"){

  }else{
    X.mod<- X.mod[,-match(excl,colnames(X.mod))]
    var.names<- var.names[match(colnames(X.mod),var.names)]
  }
  X.mod <- X.mod[,-match(names(ivs),colnames(X.mod))]
  var.names <- var.names[-match(names(ivs),var.names)]
  if(class(X.mod)=="numeric"){
    controls <- mean(X.mod)
    names(controls) <-var.names
  }else{
    controls <- apply(X.mod,2,FUN="mean")}
  design <- ivs
  design <- cbind(design,t(controls))
  return(design)}

margins.dat <- function(mod,des,alpha=.05,rounded=3,cumulate="no",pscl.data=data,num.sample=1000,prop.sample=.9){
  # Supported models include lm, glm, polr, multinom, vgam,
  # zeroinf (pscl), hurdle (pscl) and zerotrunc (countreg)
  # vgam is supported for partial proportional odds models, not models for count outcomes
  # zerotrunc is supported with bootstrapped inference, may take a while
  require(emmeans)
  if (cumulate=="no"){

    if(sum(class(mod)=="lm")>0){
      probs<-data.frame(emmeans(mod,~1,type="response",weights="proportional",at=as.list(des[1,])))
      if(nrow(des)>1){
        for(i in 2:nrow(des)){
          probsi<-data.frame(emmeans(mod,~1,type="response",weights="proportional",at=as.list(des[i,])))
          probs <- rbind(probs,probsi)}}
      probs<-probs[,-c(1,4)]
      probs[,3] <- probs[,1] - qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,2]
      probs[,4] <- probs[,1] + qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,2]
      colnames(probs)<-c("fitted","SE","ll","ul")
      marginsdat <- data.frame(des,probs)
      marginsdat<-round(marginsdat,rounded)
    }
    if(class(mod)[1] == "polr" ){
      probs<-data.frame(emmeans(mod,~dv,at=as.list(des[1,]),weights="proportional",mode="prob"))
      if(nrow(des)>1){
        for(i in 2:nrow(des)){
          probsi<-data.frame(emmeans(mod,~dv,at=as.list(des[i,]),weights="proportional",mode="prob"))
          probs <- rbind(probs,probsi)}}
      probs[,5] <- probs[,2] - qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,3]
      probs[,6] <- probs[,2] + qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,3]
      colnames(probs)[5:6]<-c("ll","ul")
      probs<-probs[,-4]
      probs[,2]<- round(probs[,2],rounded)
      probs[,3]<-round(probs[,3],rounded)
      probs[,4]<-round(probs[,4],rounded)
      probs[,5]<-round(probs[,5],rounded)
      #probs<-round(probs[,2:ncol(probs)],rounded)
      des<-des[rep(1:nrow(des), each=length(unique(probs$dv))),]
      des<-round(des,rounded)
      marginsdat <- data.frame(des,probs)}
    if(class(mod)[1] == "multinom"){
      probs<-data.frame(emmeans(mod,~dv,at=as.list(des[1,]),weights="proportional",mode="prob"))
      if(nrow(des)>1){
        for(i in 2:nrow(des)){
          probsi<-data.frame(emmeans(mod,~dv,at=as.list(des[i,]),weights="proportional",mode="prob"))
          probs <- rbind(probs,probsi)}}
      probs[,5] <- probs[,2] - qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,3]
      probs[,6] <- probs[,2] + qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,3]
      colnames(probs)[5:6]<-c("ll","ul")
      probs<-probs[,-4]
      probs[,2]<-round(probs[,2],rounded)
      probs[,3]<-round(probs[,3],rounded)
      probs[,4]<-round(probs[,4],rounded)
      probs[,5]<-round(probs[,5],rounded)
      #probs<-round(probs[,2:ncol(probs)],rounded)
      des<-des[rep(1:nrow(des), each=length(unique(probs$dv))),]
      des<-round(des,rounded)
      marginsdat <- data.frame(des,probs)}
    if(class(mod)[1] =="vglm"){
      preds<- predict(mod,newdata=des[1,],type="link",se.fit=TRUE)
      pdes <- des[rep(1,each=length(preds$fitted.values)+1),]
      pdes <- round(pdes,rounded)
      pdes <- data.frame(pdes,dv=1:(length(preds$fitted.values)+1),c(t(preds$fitted.values),NA),c(t(preds$se.fit),NA))
      colnames(pdes)[(ncol(pdes)-1):ncol(pdes)] <- c("fitted","se")
      pdes <- mutate(pdes,pr=plogis(fitted),
                     ll=plogis(fitted-qnorm(1-(alpha/2),lower.tail=TRUE)*se),
                     ul=plogis(fitted+qnorm(1-(alpha/2),lower.tail=TRUE)*se),
                     SE=(pr - ll)/qnorm(1-(alpha/2)))
      pdes <- pdes[,-match(c("fitted","se","ll","ul"),colnames(pdes))]
      pdes <- mutate(pdes,prob=NA)
      pdes$prob[1] <- pdes$pr[1]
      pdes$prob[nrow(pdes)]<- 1 - pdes$pr[nrow(pdes)-1]
      for(i in 2:(nrow(pdes)-1)) {
        pdes$prob[i] <- pdes$pr[i] - pdes$pr[i-1]
      }
      pdes <- mutate(pdes,se=NA)
      pdes$se[1] <- pdes$SE[1]
      pdes$se[nrow(pdes)]<- pdes$SE[nrow(pdes)-1]
      for(i in 2:(nrow(pdes)-1)) {
        pdes$se[i] <- sqrt(pdes$SE[i]^2 +  pdes$SE[i]^2)
      }
      pdes <- pdes[,-match(c("pr","SE"),colnames(pdes))]
      if (nrow(des)>1){for(i in 2:nrow(des)){
        predsi<- predict(mod,newdata=des[i,],type="link",se.fit=TRUE)
        pdesi <- des[rep(i,each=length(predsi$fitted.values)+1),]
        pdesi <- data.frame(pdesi,dv=1:(length(predsi$fitted.values)+1),c(t(predsi$fitted.values),NA),c(t(predsi$se.fit),NA))
        colnames(pdesi)[(ncol(pdesi)-1):ncol(pdesi)] <- c("fitted","se")
        pdesi <- mutate(pdesi,pr=plogis(fitted),
                        ll=plogis(fitted-qnorm(1-(alpha/2),lower.tail=TRUE)*se),
                        ul=plogis(fitted+qnorm(1-(alpha/2),lower.tail=TRUE)*se),
                        SE=(pr - ll)/qnorm(1-(alpha/2)))
        pdesi <- pdesi[,-match(c("fitted","se","ll","ul"),colnames(pdesi))]
        pdesi <- mutate(pdesi,prob=NA)
        pdesi$prob[1] <- pdesi$pr[1]
        pdesi$prob[nrow(pdesi)]<- 1 - pdesi$pr[nrow(pdesi)-1]
        for(i in 2:(nrow(pdesi)-1)) {
          pdesi$prob[i] <- pdesi$pr[i] - pdesi$pr[i-1]
        }
        pdesi <- mutate(pdesi,se=NA)
        pdesi$se[1] <- pdesi$SE[1]
        pdesi$se[nrow(pdesi)]<- pdesi$SE[nrow(pdesi)-1]
        for(i in 2:(nrow(pdesi)-1)) {
          pdesi$se[i] <- sqrt(pdesi$SE[i]^2 +  pdesi$SE[i]^2)
        }
        pdesi <- pdesi[,-match(c("pr","SE"),colnames(pdesi))]
        pdes<-rbind(pdes,pdesi)
      }}
      pdes <- mutate(pdes,ll=prob-qnorm(1-(alpha/2),lower.tail=TRUE)*se,
                     ul=prob+qnorm(1-(alpha/2),lower.tail=TRUE)*se)
      pdes<-mutate(pdes,prob=round(prob,rounded),
                   se=round(se,rounded),
                   ll=round(ll,rounded),
                   ul=round(ul,rounded))
      marginsdat<-pdes
    }
    if(class(mod)[1] =="zeroinfl"){
      p1<-data.frame(emmeans(mod,~1,at=as.list(des[1,]),mode="count",data=pscl.data))
      if(nrow(des)>1){
        for(i in 2:nrow(des)){
          pi<-data.frame(emmeans(mod,~1,at=as.list(des[i,]),mode="count",data=pscl.data))
          p1 <- rbind(p1,pi)}}
      p1 <- p1[,c(2,3)]
      colnames(p1)[1] <- "fitted"
      marginsdat<-cbind(des,p1)
      marginsdat<- mutate(marginsdat,ll=fitted-qnorm(1-(alpha/2),lower.tail=TRUE)*SE,
                          ul=fitted+qnorm(1-(alpha/2),lower.tail=TRUE)*SE)
      marginsdat<-round(marginsdat,rounded)
    }
    if(class(mod)[1]=="zerotrunc"){
      p1<-predict(mod,newdata=des[1,],type="response")
      if(nrow(des)>1){
        for(i in 2:nrow(des)){
          pi<-predict(mod,newdata=des[i,],type="response")
          p1<-c(p1,pi)}}
      # Boot sampling dist
      p1.model<-mod$model
      p1.dist <-matrix(NA,nr=num.sample,nc=length(p1))
      for(i in 1:num.sample){
        p1.model2 <- p1.model[sample(1:nrow(p1.model),round(prop.sample*nrow(p1.model),0),replace=TRUE),]
        p1.mod <- zerotrunc(mod$formula,data=p1.model2,dist=mod$dist)
        p1<-predict(p1.mod,newdata=des[1,],type="response")
        if(nrow(des)>1){
          for(j in 2:nrow(des)){
            p1.model2 <- p1.model[sample(1:nrow(p1.model),round(prop.sample*nrow(p1.model),0),replace=TRUE),]
            p1.mod <- zerotrunc(mod$formula,data=p1.model2,dist=mod$dist)
            pi<-predict(p1.mod,newdata=des[j,],type="response")
            p1<-c(p1,pi)}}
        p1.dist[i,]<-p1}
      p1.dist[,1]<-sort(p1.dist[,1])
      if(ncol(p1.dist)>1){for(i in 2:ncol(p1.dist)){
        p1.dist[,i]<-sort(p1.dist[,i])}}
      se <- apply(p1.dist,2,FUN="sd")
      marginsdat<-data.frame(round(des,rounded),fitted=round(p1,rounded),
                             se=round(se,rounded),
                             ll=round(p1.dist[nrow(p1.dist)*(alpha/2),],rounded),
                             ul=round(p1.dist[nrow(p1.dist)*(1-(alpha/2)),],rounded))
    }
    if(class(mod)[1]=="hurdle"){
      p1<-data.frame(emmeans(mod,~1,at=as.list(des[1,]),mode="response"))[2:3]
      if(nrow(des>1)){for(i in 2:nrow(des)){
        pi<-data.frame(emmeans(mod,~1,at=as.list(des[i,]),mode="response"))[2:3]
        p1<-rbind(p1,pi)}}
      ll <- p1$emmean - qnorm(1-(alpha/2),lower.tail=TRUE)*p1$SE
      ul <- p1$emmean + qnorm(1-(alpha/2),lower.tail=TRUE)*p1$SE
      marginsdat<-data.frame(round(des,rounded),fitted=round(p1$emmean,rounded),
                             se=round(p1$SE,rounded),
                             ll=round(ll,rounded),
                             ul=round(ul,rounded))
    }

  }else{ # cumulative probabilities
    if(class(mod)=="polr"){
      probs<-data.frame(emmeans(mod,~cut,at=as.list(des[1,]),weights="proportional",mode="cum.prob"))
      if(nrow(des)>1){
        for(i in 2:nrow(des)){
          probsi<-data.frame(emmeans(mod,~cut,at=as.list(des[i,]),weights="proportional",mode="cum.prob"))
          probs <- rbind(probs,probsi)}}
      probs[,5] <- probs[,2] - qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,3]
      probs[,6] <- probs[,2] + qnorm(1-(alpha/2),lower.tail=TRUE)*probs[,3]
      colnames(probs)[5:6]<-c("ll","ul")
      probs<-probs[,-4]
      probs[,2]<-round(probs[,2],rounded)
      probs[,3]<-round(probs[,3],rounded)
      probs[,4]<-round(probs[,4],rounded)
      probs[,5]<-round(probs[,5],rounded)
      des<-des[rep(1:nrow(des), each=length(table(mod$model[,1]))-1),]
      des<-round(des,rounded)
      marginsdat <- data.frame(des,probs)}
    if(class(mod)=="vglm"){
      preds<- predict(mod,newdata=des[1,],type="link",se.fit=TRUE)
      pdes <- des[rep(1,each=length(preds$fitted.values)),]
      pdes <- round(pdes,rounded)
      pdes <- data.frame(pdes,colnames(preds$fitted.values),t(preds$fitted.values),t(preds$se.fit))
      colnames(pdes)[(ncol(pdes)-2):ncol(pdes)] <- c("dv","fitted","se")
      if(nrow(des)>1){for(i in 2:nrow(des)){
        predsi<- predict(mod,newdata=des[i,],type="link",se.fit=TRUE)
        pdesi <- des[rep(i,each=length(predsi$fitted.values)),]
        pdesi <- round(pdesi,rounded)
        pdesi <- data.frame(pdesi,colnames(preds$fitted.values),t(predsi$fitted.values),t(predsi$se.fit))
        colnames(pdesi)[(ncol(pdesi)-2):ncol(pdesi)] <- c("dv","fitted","se")
        pdes <- rbind(pdes,pdesi)}}
      pdes <- mutate(pdes,prob=plogis(fitted),
                     ll=plogis(fitted-qnorm(1-(alpha/2),lower.tail=TRUE)*se),
                     ul=plogis(fitted+qnorm(1-(alpha/2),lower.tail=TRUE)*se),
                     SE=(prob - ll)/qnorm(1-(alpha/2)))
      pdes <- pdes[,-match(c("fitted","se"),colnames(pdes))]
      pdes<-mutate(pdes,prob=round(prob,rounded),
                   SE=round(SE,rounded),
                   ll=round(ll,rounded),
                   ul=round(ul,rounded))
      marginsdat<-pdes}
  } # end cumulative probabilities
  colnames(marginsdat)[match("SE",colnames(marginsdat))] <- "se"
  return(marginsdat)
}

margins.dat.clogit <- function(mod,design.matrix,run.boot="no",num.sample=1000,prop.sample=.9,alpha=.05){
  require(tidyverse)
  des <- mutate(design.matrix,lp=exp(as.matrix(design.matrix) %*% coef(mod)),
                probs=lp/sum(lp))
  out<-des
  boot.dist <- matrix(NA,nr=num.sample,nc=nrow(des))
  if(run.boot=="yes"){
    for(i in 1:num.sample){
      mod2 <- mod$model[sample(1:nrow(mod$model),round(prop.sample*nrow(mod$model),0),replace=TRUE),]
      m.1 <- clogistic(mod$formula, strata=`(strata)`, data=mod2)
      design2 <- mutate(des,lp=exp(as.matrix(design) %*% coef(m.1)),
                        probs=lp/sum(lp))
      boot.dist[i,] <-design2[,ncol(design2)]
    }
    boot.dist[,1] <- sort(boot.dist[,1])
    if(ncol(boot.dist)>1){
      for(i in 2:ncol(boot.dist)){
        boot.dist[,i] <- sort(boot.dist[,i])
      }
    }
    lower.limit <- boot.dist[nrow(boot.dist)*(alpha/2),]
    upper.limit <- boot.dist[nrow(boot.dist)*(1-(alpha/2)),]
    des<-mutate(des,ll=lower.limit,ul=upper.limit)
    out<-list(des=des,boot.dist=boot.dist)
  }
  return(out)
}


# Take a Poisson model object of the count process.
count.fit<-function(m1,y.range,rounded=3,use.color="yes"){
  require(MASS)
  require(ggpubr)
  require(pscl)
  outcome<-m1$model[,1]
  form <- as.character(m1$formula)
  form.2<-as.formula(paste(form[2],form[1],form[3],"|",form[3]))
  m1<-glm(m1$formula ,family="poisson",data=m1$data )
  m2<-glm.nb(m1$formula ,data=m1$data )
  m3 <- zeroinfl(form.2,data=m1$data)
  m4 <- zeroinfl(form.2,data=m1$data, dist="negbin")
  bics<-c(BIC(m1),BIC(m2),BIC(m3),BIC(m4))
  aics<-c(AIC(m1),AIC(m2),AIC(m3),AIC(m4))
  out<-matrix(0,nr=2*(ncol(m1$model)),nc=16)
  out[,1]<-c(coef(m1),rep(0,length(coef(m1))))
  out[,5]<-c(coef(m2),rep(0,length(coef(m2))))
  out[,9]<-c(coef(m3))
  out[,13]<-c(coef(m4))
  out[,2]<-c(sqrt(diag(vcov((m1)))),rep(0,length(coef(m1))))
  out[,6]<-c(sqrt(diag(vcov((m2)))),rep(0,length(coef(m2))))
  out[,10]<-c(sqrt(diag(vcov((m3)))))
  out[,14]<-c(sqrt(diag(vcov((m4)))))
  out[,3]<-out[,1]/out[,2]
  out[,7]<-out[,5]/out[,6]
  out[,11]<-out[,9]/out[,10]
  out[,15]<-out[,13]/out[,14]
  out[,4]<-dnorm(out[,3])
  out[,8]<-dnorm(out[,7])
  out[,12]<-dnorm(out[,11])
  out[,16]<-dnorm(out[,15])
  out<-round(out,rounded)
  rownames(out)<-names(coef(m3))
  colnames(out)<-c("Pcoef","Pse","Pz","Ppval","NBcoef","NBse","NBz","NBpval","ZIPcoef","ZIPse","ZIPz","ZIPpval","ZNBcoef","ZNBse","ZNBz","ZNBpval")
  # Observed probabilities
  out2<-outcome[which(outcome<(max(y.range)+1))]
  obs<-0:max(y.range)
  for(i in 1:length(y.range)){
    obs[i]<-sum((out2==y.range[i]))/length(out2)}

  poi<-round(dpois(y.range,lambda=mean(predict(m1,newdata=,type="response"))),rounded)
  nb<-round(dnbinom(y.range,m2$theta,mu=mean(predict(m2,newdata=,type="response"))),rounded)
  zpoi<-round(apply(predprob(m3),2,mean)[1:length(y.range)],rounded)
  znb<-round(apply(predprob(m4),2,mean)[1:length(y.range)],rounded)
  poi<-obs-poi
  nb<-obs-nb
  zpoi<-obs-zpoi
  znb<-obs-znb
  type<-c(rep("Pois",length(y.range)),rep("NB",length(y.range)),rep("ZIP",length(y.range)),rep("ZNB",length(y.range)))
  dat<-data.frame(c(poi,nb,zpoi,znb),rep(y.range,4),type)

  if(use.color=="yes"){
    out.plot<-ggplot(dat, aes(x=dat[,2], y=dat[,1]),group=dat[,3],linetype=dat[,3]) +
      geom_hline(yintercept=0,color="gray80") +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x="Outcome",y="Observed - Predicted",color="") +
      scale_x_continuous(breaks=seq(min(y.range),max(y.range),round(length(y.range)/4)))  +
      geom_line(aes(linetype=type,color=type)) + labs(linetype="",color="") +
      scale_shape_manual(name="",values=c(22,21,23,24)) +
      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))}else{
        out.plot<-ggplot(dat, aes(x=dat[,2], y=dat[,1]),group=dat[,3],linetype=dat[,3]) +
          geom_hline(yintercept=0,color="gray80") +
          theme_bw() +
          theme(legend.position="bottom") +
          xlab("Outcome") +
          ylab("Observed - Predicted") +
          scale_x_continuous(breaks=seq(min(y.range),max(y.range),round(length(y.range)/4)))  +
          geom_line(aes(linetype=type)) + labs(linetype="") +
          scale_shape_manual(name="",values=c(22,21,23,24)) +
          theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

      }

  models <- out[,c(1,2,5,6,9,10,13,14)]
  pdat <- data.frame(var=c("Intercept",names(m1$coef)[2:length(names(m1$coefficients))]),
                     coef=as.numeric(c(models[,1],models[,3],models[,5],models[,7])),
                     se=as.numeric(c(models[,2],models[,4],models[,6],models[,8])),
                     model=rep(c("Poisson","Neg Binom","ZIP","ZNB"),each=nrow(models)))

  pdat <- mutate(pdat,max=coef+1.96*se,min=coef-1.96*se)

  pdat$var <- paste(rep(rep(c("c","z"),each=length(m1$coefficients)),4),pdat$var)
  pdat1<-pdat[1:(nrow(pdat)/4),]
  pdat1$var <- factor(pdat1$var,levels=paste(pdat1$var))
  pdat1[((nrow(pdat1)/2)+1):nrow(pdat1),2:3]<-NA
  p1<-ggplot(pdat1,aes(y=var,x=coef,xmin=coef-1.96*se,xmax=coef+1.96*se)) +
    geom_vline(xintercept=0,color="gray80") +
    geom_pointrange() +
    theme_classic() + labs(x="",y="",title="Poisson")  + scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))
  pdat2<-pdat[((nrow(pdat)/4)+1):((nrow(pdat)/2)),]
  pdat2$var <- factor(pdat2$var,levels=paste(pdat2$var))
  pdat2[((nrow(pdat2)/2)+1):nrow(pdat2),2:3]<-NA
  p2<-ggplot(pdat2,aes(y=var,x=coef,xmin=coef-1.96*se,xmax=coef+1.96*se)) +
    geom_vline(xintercept=0,color="gray80") +
    geom_pointrange() +
    theme_classic() + labs(x="",y="",title="Neg Binom")  + scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))
  pdat3<-pdat[((nrow(pdat)/2)+1):(((nrow(pdat)/4))*3),]
  pdat3$var <- factor(pdat3$var,levels=paste(pdat3$var))
  p3<-ggplot(pdat3,aes(y=var,x=coef,xmin=coef-1.96*se,xmax=coef+1.96*se)) +
    geom_vline(xintercept=0,color="gray80") +
    geom_pointrange() +
    theme_classic() + labs(x="",y="",title="ZIP")  + scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))
  pdat4<-pdat[(1+(((nrow(pdat)/4))*3)):nrow(pdat),]
  pdat4$var <- factor(pdat4$var,levels=paste(pdat4$var))
  p4<-ggplot(pdat4,aes(y=var,x=coef,xmin=coef-1.96*se,xmax=coef+1.96*se)) +
    geom_vline(xintercept=0,color="gray80") +
    geom_pointrange() +
    theme_classic() + labs(x="",y="",title="ZINB")  + scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))
  coef.pic<-ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)

  information.criterion<-rbind(bics,aics)
  rownames(information.criterion)<-c("BIC","AIC")
  colnames(information.criterion)<-c("Poisson","Neg Binom","ZIP","ZNB")

  vuong(m1,m3)
  vuong(m2,m4)

  outp<-list(ic=information.criterion,pic=out.plot,models=out,models.pic=coef.pic)
  return(outp)}




first.diff.fitted  <- function(mod,design.matrix,compare,alpha=.05,rounded=3,bootstrap="no",num.sample=1000,prop.sample=.9,data){
  if (bootstrap=="no"){
    require(marginaleffects)
    if(sum(class(mod)=="lm")>0){
      f1 <- function(x) predict(x, type = "response", newdata = design.matrix[compare[1],]) - predict(x, type = "response", newdata = design.matrix[compare[2],])
    }else{
      f1 <- function(x) predict(x, type = "probs", newdata = design.matrix[compare[1],]) - predict(x, type = "probs", newdata = design.matrix[compare[2],])
    }
    out <- deltamethod(mod, FUN = f1)
    colnames(out)[c(2,6:7)]<-c("est","ll","ul")
    out[,6] <- out[,2]-qnorm(1-(alpha/2),lower.tail=TRUE)*out[,3]
    out[,7] <- out[,2]+qnorm(1-(alpha/2),lower.tail=TRUE)*out[,3]

    out[,2:7]<-round(out[,2:7],rounded)
  }
  if (bootstrap=="yes"){
    m8 <- as.character(mod$call)

    if(m8[1]=="lm") {
      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- lm(m8[2],data=dat2)
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- lm(m8[2],data=dat2)
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class="OLS"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,ci.95=ci.95,model.class=model.class)
    }

    if(m8[1]=="glm" & m8[3]=="binomial") {

      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- glm(m8[2],data=dat2,family="binomial")
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- glm(m8[2],data=dat2,family="binomial")
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class="logit"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,ci.95=ci.95,model.class=model.class)
    }
    if (m8[1]=="glm" & m8[3]=="poisson") {
      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- glm(m8[2],data=dat2,family="poisson")
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- glm(m8[2],data=dat2,family="poisson")
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class <- "Poisson"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)

    }
    if (m8[1]=="glm.nb") {
      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- glm.nb(m8[2],data=dat2)
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- glm.nb(m8[2],data=dat2)
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class<-"negative binomial"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)
    }
    if (m8[1]=="polr")  {

      obs.diff <- as.numeric(predict(mod, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod, type = "probs", newdata = design.matrix[compare[2],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2<-polr(m8[2],data=dat2, Hess=TRUE)
      d1 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- polr(m8[2],data=dat2, Hess=TRUE)
        d2 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],]))
        d1 <- rbind(d1,d2)}
      boot.dist <- d1
      boot.dist[,1] <- sort(boot.dist[,1])
      for(i in 2:ncol(boot.dist)){
        boot.dist[,i] <- sort(boot.dist[,i])
      }
      mean.boot.dist <- apply(boot.dist,2,FUN="mean")
      sd.boot.dist <- apply(boot.dist,2,FUN="sd")
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),1]
      for(i in 2:ncol(boot.dist)){
        ci.952 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),i]
        ci.95<- rbind(ci.95,ci.952)}
      colnames(ci.95) <- c("Lower Limit","Upper Limit")
      model.class <- "ordered logit"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)
    }
    if (m8[1]=="multinom"){
      obs.diff <- as.numeric(predict(mod, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod, type = "probs", newdata = design.matrix[compare[2],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2<-multinom(m8[2],data=dat2)

      d1 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2<-multinom(m8[2],data=dat2)
        d2 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],]))
        d1 <- rbind(d1,d2)}
      boot.dist <- d1
      boot.dist[,1] <- sort(boot.dist[,1])
      for(i in 2:ncol(boot.dist)){
        boot.dist[,i] <- sort(boot.dist[,i])
      }
      mean.boot.dist <- apply(boot.dist,2,FUN="mean")
      sd.boot.dist <- apply(boot.dist,2,FUN="sd")
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),1]
      for(i in 2:ncol(boot.dist)){
        ci.952 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),i]
        ci.95<- rbind(ci.95,ci.952)}
      colnames(ci.95) <- c("Lower Limit","Upper Limit")
      model.class <- "multinomial"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)

    }


  }
  return(out)  }
second.diff.fitted  <- function(mod,design.matrix,compare,alpha=.05,rounded=3,bootstrap="no",num.sample=1000,prop.sample=.9,data){

  if(bootstrap=="no"){
    require(marginaleffects)
    if(sum(class(mod)=="lm")>0){
      f1 <- function(x) (predict(x, type = "response", newdata = design.matrix[compare[1],]) - predict(x, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(x, type = "response", newdata = design.matrix[compare[3],]) - predict(x, type = "response", newdata = design.matrix[compare[4],]))
    }else{
      f1 <- function(x) (predict(x, type = "probs", newdata = design.matrix[compare[1],]) - predict(x, type = "probs", newdata = design.matrix[compare[2],])) -
        (predict(x, type = "probs", newdata = design.matrix[compare[3],]) - predict(x, type = "probs", newdata = design.matrix[compare[4],]))
    }
    out <- deltamethod(mod, FUN = f1)
    colnames(out)[c(2,6:7)]<-c("est","ll","ul")
    out[,6] <- out[,2]-qnorm(1-(alpha/2),lower.tail=TRUE)*out[,3]
    out[,7] <- out[,2]+qnorm(1-(alpha/2),lower.tail=TRUE)*out[,3]
    out[,2:7]<-round(out[,2:7],rounded)
  }
  if(bootstrap=="yes"){
    m8 <- as.character(mod$call)


    if(m8[1]=="lm") {
      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod, type = "response", newdata = design.matrix[compare[3],]) - predict(mod, type = "response", newdata = design.matrix[compare[4],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- lm(m8[2],data=dat2)
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- lm(m8[2],data=dat2)
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
          (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class="OLS"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,ci.95=ci.95,model.class=model.class)
    }

    if(m8[1]=="glm" & m8[3]=="binomial") {

      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod, type = "response", newdata = design.matrix[compare[3],]) - predict(mod, type = "response", newdata = design.matrix[compare[4],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- glm(m8[2],data=dat2,family="binomial")
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- glm(m8[2],data=dat2,family="binomial")
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
          (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class="logit"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,ci.95=ci.95,model.class=model.class)
    }
    if (m8[1]=="glm" & m8[3]=="poisson") {
      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod, type = "response", newdata = design.matrix[compare[3],]) - predict(mod, type = "response", newdata = design.matrix[compare[4],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- glm(m8[2],data=dat2,family="poisson")
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- glm(m8[2],data=dat2,family="poisson")
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
          (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class <- "Poisson"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)

    }
    if (m8[1]=="glm.nb") {
      obs.diff <- as.numeric(predict(mod, type = "response", newdata = design.matrix[compare[1],]) - predict(mod, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod, type = "response", newdata = design.matrix[compare[3],]) - predict(mod, type = "response", newdata = design.matrix[compare[4],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2 <- glm.nb(m8[2],data=dat2)
      d1<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
        (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- glm.nb(m8[2],data=dat2)
        d2<- as.numeric(predict(mod2, type = "response", newdata = design.matrix[compare[1],]) - predict(mod2, type = "response", newdata = design.matrix[compare[2],])) -
          (predict(mod2, type = "response", newdata = design.matrix[compare[3],]) - predict(mod2, type = "response", newdata = design.matrix[compare[4],]))
        d1 <- c(d1,d2)}

      boot.dist <- sort(d1)
      mean.boot.dist <- mean(boot.dist)
      sd.boot.dist <- sd(boot.dist)
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample)]
      names(ci.95) <- c("Lower Limit","Upper Limit")
      model.class<-"negative binomial"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)
    }
    if (m8[1]=="polr")  {

      obs.diff <- as.numeric(predict(mod, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod, type = "probs", newdata = design.matrix[compare[2],])) -
        (predict(mod, type = "probs", newdata = design.matrix[compare[3],]) - predict(mod, type = "probs", newdata = design.matrix[compare[4],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2<-polr(m8[2],data=dat2, Hess=TRUE)
      d1 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],])) -
        (predict(mod2, type = "probs", newdata = design.matrix[compare[3],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[4],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2 <- polr(m8[2],data=dat2, Hess=TRUE)
        d2 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],])) -
          (predict(mod2, type = "probs", newdata = design.matrix[compare[3],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[4],]))
        d1 <- rbind(d1,d2)}
      boot.dist <- d1
      boot.dist[,1] <- sort(boot.dist[,1])
      for(i in 2:ncol(boot.dist)){
        boot.dist[,i] <- sort(boot.dist[,i])
      }
      mean.boot.dist <- apply(boot.dist,2,FUN="mean")
      sd.boot.dist <- apply(boot.dist,2,FUN="sd")
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),1]
      for(i in 2:ncol(boot.dist)){
        ci.952 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),i]
        ci.95<- rbind(ci.95,ci.952)}
      colnames(ci.95) <- c("Lower Limit","Upper Limit")
      model.class <- "ordered logit"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)
    }
    if (m8[1]=="multinom"){
      obs.diff <- as.numeric(predict(mod, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod, type = "probs", newdata = design.matrix[compare[2],])) -
        (predict(mod, type = "probs", newdata = design.matrix[compare[3],]) - predict(mod, type = "probs", newdata = design.matrix[compare[4],]))
      dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
      mod2<-multinom(m8[2],data=dat2)

      d1 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],])) -
        (predict(mod2, type = "probs", newdata = design.matrix[compare[3],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[4],]))
      for (i in 2:num.sample){
        dat2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
        mod2<-multinom(m8[2],data=dat2)
        d2 <- as.numeric(predict(mod2, type = "probs", newdata = design.matrix[compare[1],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[2],])) -
          (predict(mod2, type = "probs", newdata = design.matrix[compare[3],]) - predict(mod2, type = "probs", newdata = design.matrix[compare[4],]))
        d1 <- rbind(d1,d2)}
      boot.dist <- d1
      boot.dist[,1] <- sort(boot.dist[,1])
      for(i in 2:ncol(boot.dist)){
        boot.dist[,i] <- sort(boot.dist[,i])
      }
      mean.boot.dist <- apply(boot.dist,2,FUN="mean")
      sd.boot.dist <- apply(boot.dist,2,FUN="sd")
      ci.95 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),1]
      for(i in 2:ncol(boot.dist)){
        ci.952 <- boot.dist[c((alpha/2)*num.sample,(1-alpha/2)*num.sample),i]
        ci.95<- rbind(ci.95,ci.952)}
      colnames(ci.95) <- c("Lower Limit","Upper Limit")
      model.class <- "multinomial"
      out<-list(obs.diff=obs.diff,boot.dist=boot.dist,mean.boot.dist=mean.boot.dist,sd.boot.dist=sd.boot.dist,
                ci.95=ci.95,model.class=model.class)

    }

  }
  return(out)}




