
list.coef<-function(model,rounded=3,alpha=.05){

  if(is(model,"glmerMod") | is(model,"lmerMod")){
    out<-matrix(0,nrow=length(nlme::fixef(model)),ncol=10)
    }else{out<-matrix(0,nrow=length(coef(model)),ncol=10)}

  if(is(model,"multinom")){out[,1]<-t(coef(model))
  }else if(is(model,"lme")){
    out[,1]<- nlme::fixef(model)
  }else if(is(model,"glmerMod") | is(model,"lmerMod")){
      out[,1]<- nlme::fixef(model)}else{out[,1]<-coef(model)}

  if(is(model,"lme")){vcov1 <- vcov(model)}else if(is(model,"glmerMod") | is(model,"lmerMod")){vcov1 <- vcov(model)}else{vcov1<-vcov(model)[1:length(coef(model)),1:length(coef(model))]}
  out[,2]<-sqrt(diag(vcov1))
  out[,3]<-out[,1]/out[,2]
  out[,4] <- out[,1]-qnorm(1-(alpha/2),lower.tail=TRUE)*out[,2]
  out[,5] <- out[,1]+qnorm(1-(alpha/2),lower.tail=TRUE)*out[,2]
  out[,6]<-dnorm(out[,3])
  if(is(model,"multinom")){ out[,7]<-exp(t(coef(model)))}else if(is(model,"lme")){out[,7]<-exp(nlme::fixef(model))}else if(is(model,"glmerMod") | is(model,"lmerMod")){out[,7]<-exp(nlme::fixef(model))}else{out[,7]<-exp(coef(model))}
  out[,8]<-exp(out[,4])
  out[,9]<-exp(out[,5])
  out[,10]<-100*(exp(out[,1])-1)
  out<-round(out,rounded)
  colnames(out)<-c("b","SE","z","ll","ul","p.val","exp.b","ll.exp.b","ul.exp.b","percent")
  out <- data.frame(out,CI=paste(100*(1-alpha),"%"))
  if(is(model,"multinom")){
    cn<-colnames(coef(model))
    rn<-rownames(coef(model))
    names <- paste(rn[1],cn)
    for(i in 2:length(rn)){
      namesi <- paste(rn[i],cn)
      names <- c(names,namesi)}
    out <- data.frame(variables=names,out)
  }else if(is(model,"glmerMod") | is(model,"lmerMod")){out<-data.frame(variables=names(nlme::fixef(model)),out)}else{
    out <- data.frame(variables=names(coef(model)),out)}
  return(out)}



margins.des<-function (mod, ivs, excl = "nonE",data) {
  if(is(mod,"lm") | is(mod,"glm") | is(mod,"polr") | is(mod,"multinom") | is(mod,"vglm")  |
    is(mod,"negbin") | is(mod,"zeroinfl") | is(mod,"hurdle") | is(mod,"glmerMod")  |
    is(mod,"clmm")  | is(mod,"lme") | is(mod,"lmerModLmerTest") | is(mod,"lmerMod")){


    if(is(mod,"nnet")){
      c1<-as.character(mod$call)
      m.polr<-MASS::polr(c1[2],data=data)
      mod<-m.polr}

    if(is(mod,"lme")){
      X5 <- mod$data # pull out all of the data
      var.names <- names(mod$fixDF$terms) # get the names of the variables that were used
      var.names <- var.names[-1] # remove the intercept
      X2 <- X5[,match(var.names,colnames(X5))] # from the full data, keep only the variables in the model
      X.mod <- X2 # rename
      if (excl == "nonE") {}else {
        X.mod <- X2[, -match(excl, colnames(X2))]
        var.names <- var.names[match(colnames(X.mod), var.names)]} # If you want to remove any factor-coded variables

      X.mod <- X.mod[, -match(names(ivs), colnames(X.mod))] # Exclude "ivs" because those are not set to their means
      var.names <- var.names[-match(names(ivs), var.names)] # Exclude "ivs" because those are not set to their means
      if (is(X.mod,"numeric")) {
        controls <- mean(X.mod)
        names(controls) <- var.names
      }else {
        controls <- apply(X.mod, 2, FUN = "mean")} # Set covariates to their means
    }else if(is(mod,"glmerMod")){
      dv <- as.character(formula(mod))[[2]]
      var.names <- names(nlme::fixef(mod)) #names of variables in the model
      var.names <- var.names[-1] #remove intercept
      # Seems like I need to add a variable to the function that defines the DV (so I can adjust for missing data)
      data <- data[, match(c(dv,var.names), colnames(data))]
      X.mod <- na.omit(data) # remove missing responses
      X.mod <- X.mod[,-1] # Now, remove the dv from the predictors
      if (excl == "nonE") {
      }else {
        X.mod <- X.mod[, -match(excl, colnames(X.mod))]
        var.names <- var.names[match(colnames(X.mod), var.names)]
      }
      X.mod <- X.mod[, -match(names(ivs), colnames(X.mod))] # remove variables in the ivs statement
      var.names <- var.names[-match(names(ivs), var.names)]
      if (is(X.mod,"numeric")) {
        controls <- mean(X.mod)
        names(controls) <- var.names
      }else {
        controls <- apply(X.mod, 2, FUN = "mean")
      }
    }else if(is(mod,"lmerModLmerTest") | is(mod,"lmerMod")){
      dv <- as.character(formula(mod))[[2]]
      var.names <- names(nlme::fixef(mod)) #names of variables in the model
      var.names <- var.names[-1] #remove intercept
      # Seems like I need to add a variable to the function that defines the DV (so I can adjust for missing data)
      data <- data[, match(c(dv,var.names), colnames(data))]
      X.mod <- na.omit(data) # remove missing responses
      X.mod <- X.mod[,-1] # Now, remove the dv from the predictors
      if (excl == "nonE") {
      }else {
        X.mod <- X.mod[, -match(excl, colnames(X.mod))]
        var.names <- var.names[match(colnames(X.mod), var.names)]
      }
      X.mod <- X.mod[, -match(names(ivs), colnames(X.mod))] # remove variables in the ivs statement
      var.names <- var.names[-match(names(ivs), var.names)]
      if (is(X.mod,"numeric")) {
        controls <- mean(X.mod)
        names(controls) <- var.names
      }else {
        controls <- apply(X.mod, 2, FUN = "mean")
      }
    }else if(is(mod,"clmm")){
      X <- mod$model
      X <- X[,-c(1,ncol(X))]
      controls <- apply(X, 2, FUN = "mean")
      controls <- controls[-match(names(ivs), colnames(X))]
    }
    else{
      X.mod <- mod$model[, -1]
      var.names <- colnames(X.mod)
      if (excl == "nonE") {
      }
      else {
        X.mod <- X.mod[, -match(excl, colnames(X.mod))]
        var.names <- var.names[match(colnames(X.mod), var.names)]
      }
      X.mod <- X.mod[, -match(names(ivs), colnames(X.mod))]
      var.names <- var.names[-match(names(ivs), var.names)]
      if (is(X.mod,"numeric")) {
        controls <- mean(X.mod)
        names(controls) <- var.names
      }
      else {
        controls <- apply(X.mod, 2, FUN = "mean")
      }}

    design <- ivs
    design <- cbind(design, t(controls))
    return(design)} else {message("Model type not supported.")}
}

margins.dat <- function (mod, des, alpha = 0.05, rounded = 3, cumulate = "no",
                         pscl.data = data, num.sample = 1000, prop.sample = 0.9, seed = 1234) {
  pr <- NULL; prob <- NULL; SE <- NULL; fitted <- NULL ; se <- NULL
  if(is(mod,"lm") | is(mod,"glm") | is(mod,"polr") | is(mod,"multinom") | is(mod,"vglm")  |
     is(mod,"negbin")  | is(mod,"zeroinfl") |  is(mod,"hurdle") | is(mod,"glmerMod")  |
     is(mod,"clmm")  | is(mod,"lme")| is(mod,"lmerModLmerTest") | is(mod,"lmerMod")){


    if (cumulate == "no") {
      if (sum(class(mod) == "lm") > 0) { # if I have to change this line, I'll need to differentiate various forms of glms
        probs <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",
                                                              weights = "proportional", at = as.list(des[1,
                                                              ]))))
        if (nrow(des > 1)) {
          for (i in 2:nrow(des)) {
            probsi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",
                                                                   weights = "proportional", at = as.list(des[i,
                                                                   ]))))
            probs <- rbind(probs, probsi)
          }
        }
        probs <- probs[, -c(1, 4)]
        probs[, 3] <- probs[, 1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 2]
        probs[, 4] <- probs[, 1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 2]
        colnames(probs) <- c("fitted", "SE", "ll", "ul")
        marginsdat <- data.frame(des, probs)
        marginsdat <- round(marginsdat, rounded)
      }
      if (is(mod,"polr")) {
        probs <- suppressMessages(data.frame(emmeans::emmeans(mod, specs=colnames(mod$model)[1], at = as.list(des[1,
        ]), weights = "proportional", mode = "prob")))
        if (nrow(des) > 1) {
          for (i in 2:nrow(des)) {
            probsi <- suppressMessages(data.frame(emmeans::emmeans(mod, specs=colnames(mod$model)[1], at = as.list(des[i,
            ]), weights = "proportional", mode = "prob")))
            probs <- rbind(probs, probsi)
          }
        }
        probs[, 5] <- probs[, 2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 3]
        probs[, 6] <- probs[, 2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 3]
        colnames(probs)[5:6] <- c("ll", "ul")
        probs <- probs[, -4]
        probs[, 2] <- round(probs[, 2], rounded)
        probs[, 3] <- round(probs[, 3], rounded)
        probs[, 4] <- round(probs[, 4], rounded)
        probs[, 5] <- round(probs[, 5], rounded)
        des <- des[rep(1:nrow(des), each = length(table(mod$model[,1]))),
        ]
        des <- round(des, rounded)
        marginsdat <- data.frame(des, probs)
      }
      if (is(mod,"multinom")) {
        cl <- as.character(mod$call)[2]
        dv<-strsplit(cl," ")[[1]][1]
        probs <- suppressMessages(data.frame(emmeans::emmeans(mod, specs=dv, at = as.list(des[1,
        ]), weights = "proportional", mode = "prob")))
        if (nrow(des) > 1) {
          for (i in 2:nrow(des)) {
            probsi <- suppressMessages(data.frame(emmeans::emmeans(mod, specs=dv, at = as.list(des[i,
            ]), weights = "proportional", mode = "prob")))
            probs <- rbind(probs, probsi)
          }
        }
        probs[, 5] <- probs[, 2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 3]
        probs[, 6] <- probs[, 2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 3]
        colnames(probs)[5:6] <- c("ll", "ul")
        probs <- probs[, -4]
        probs[, 2] <- round(probs[, 2], rounded)
        probs[, 3] <- round(probs[, 3], rounded)
        probs[, 4] <- round(probs[, 4], rounded)
        probs[, 5] <- round(probs[, 5], rounded)
        des <- des[rep(1:nrow(des), each = length(mod$lev)),
        ]
        des <- round(des, rounded)
        marginsdat <- data.frame(des, probs)
      }



      if (is(mod,"zeroinfl")) {
        p1 <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, at = as.list(des[1,
        ]), mode = "count", data = pscl.data)))
        if (nrow(des) > 1) {
          for (i in 2:nrow(des)) {
            pi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, at = as.list(des[i,
            ]), mode = "count", data = pscl.data)))
            p1 <- rbind(p1, pi)
          }
        }
        p1 <- p1[, c(2, 3)]
        colnames(p1)[1] <- "fitted"
        marginsdat <- cbind(des, p1)
        marginsdat <- dplyr::mutate(marginsdat, ll = fitted - qnorm(1 -
                                                                      (alpha/2), lower.tail = TRUE) * SE, ul = fitted +
                                      qnorm(1 - (alpha/2), lower.tail = TRUE) * SE)
        marginsdat <- round(marginsdat, rounded)
      }

      if (is(mod,"hurdle")) {
        p1 <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, at = as.list(des[1,
        ]), mode = "response"))[2:3])
        if (nrow(des > 1)) {
          for (i in 2:nrow(des)) {
            pi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, at = as.list(des[i,
            ]), mode = "response"))[2:3])
            p1 <- rbind(p1, pi)
          }
        }
        ll <- p1$emmean - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        ul <- p1$emmean + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        marginsdat <- data.frame(round(des, rounded), fitted = round(p1$emmean,
                                                                     rounded), se = round(p1$SE, rounded), ll = round(ll,
                                                                                                                      rounded), ul = round(ul, rounded))
      }

      if (is(mod,"lme")) {
        p1 <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",
                                                           weights = "proportional", at = as.list(des[1,])))[2:3])
        if (nrow(des > 1)) {
          for (i in 2:nrow(des)) {
            pi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",
                                                               weights = "proportional", at = as.list(des[i,])))[2:3])
            p1 <- rbind(p1, pi)
          }
        }
        ll <- p1$emmean - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        ul <- p1$emmean + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        marginsdat <- data.frame(round(des, rounded), fitted = round(p1$emmean,
                                                                     rounded), se = round(p1$SE, rounded), ll = round(ll,
                                                                                                                      rounded), ul = round(ul, rounded))

      }

      if (is(mod,"glmerMod")) {
        p1 <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",weights = "proportional", at = as.list(des[1,]))))
        if (nrow(des > 1)) {
          for (i in 2:nrow(des)) {
            pi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",weights = "proportional", at = as.list(des[i,]))))
            p1 <- rbind(p1, pi)
          }
        }
        ll <- p1[,2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        ul <- p1[,2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        marginsdat <- data.frame(round(des, rounded), fitted = round(p1[,2],rounded), se = round(p1[,3], rounded), ll = round(ll,rounded), ul = round(ul, rounded))
      }
      if (is(mod,"lmerModLmerTest") | is(mod,"lmerMod")){
        p1 <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",weights = "proportional", at = as.list(des[1,]))))
        if (nrow(des > 1)) {
          for (i in 2:nrow(des)) {
            pi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~1, type = "response",weights = "proportional", at = as.list(des[i,]))))
            p1 <- rbind(p1, pi)
          }
        }
        ll <- p1[,2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        ul <- p1[,2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          p1$SE
        marginsdat <- data.frame(round(des, rounded), fitted = round(p1[,2],rounded), se = round(p1[,3], rounded), ll = round(ll,rounded), ul = round(ul, rounded))
      }
      if (is(mod,"clmm")) {
        out <- suppressMessages(data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(des[1,]))))
        if(nrow(des) >1){for(i in 2:nrow(des)){
          out2 <- suppressMessages(data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(des[i,]))))
          out <- rbind(out,out2)}}
        ll <- out[,2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          out$SE
        ul <- out[,2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          out$SE
        marginsdat <- data.frame(des,dv.level = out[,1], fitted = round(out[,2],rounded), se = round(out[,3], rounded), ll = round(ll,rounded), ul = round(ul, rounded))
      }
    }
    else {
      if (is(mod,"polr")) {
        probs <- suppressMessages(data.frame(emmeans::emmeans(mod, ~cut, at = as.list(des[1,
        ]), weights = "proportional", mode = "cum.prob")))
        if (nrow(des) > 1) {
          for (i in 2:nrow(des)) {
            probsi <- suppressMessages(data.frame(emmeans::emmeans(mod, ~cut, at = as.list(des[i,
            ]), weights = "proportional", mode = "cum.prob")))
            probs <- rbind(probs, probsi)
          }
        }
        probs[, 5] <- probs[, 2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 3]
        probs[, 6] <- probs[, 2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          probs[, 3]
        colnames(probs)[5:6] <- c("ll", "ul")
        probs <- probs[, -4]
        probs[, 2] <- round(probs[, 2], rounded)
        probs[, 3] <- round(probs[, 3], rounded)
        probs[, 4] <- round(probs[, 4], rounded)
        probs[, 5] <- round(probs[, 5], rounded)
        des <- des[rep(1:nrow(des), each = length(table(mod$model[,
                                                                  1])) - 1), ]
        des <- round(des, rounded)
        marginsdat <- data.frame(des, probs)
      }
      if (is(mod,"clmm")) {
        out <- suppressMessages(data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(des[1,]))))
        if(nrow(des) >1){for(i in 2:nrow(des)){
          out2 <- suppressMessages(data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(des[i,]))))
          out <- rbind(out,out2)}}
        ll <- out[,2] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
          out$SE
        ul <- out[,2] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
          out$SE
        marginsdat <- data.frame(round(des, rounded),cut=out[,1], fitted = round(out[,2],rounded), se = round(out[,3], rounded), ll = round(ll,rounded), ul = round(ul, rounded))
      }
    }
    colnames(marginsdat)[match("SE", colnames(marginsdat))] <- "se"
    return(marginsdat)
  } else {message("Model type not supported.")}}

# Take a Poisson model object of the count process.
count.fit<-function(m1,y.range,rounded=3,use.color="yes"){
  se <- NULL; a <- NULL; b <- NULL
  if(is(m1,"glm") & m1$family[1]$family=="poisson" ){
    outcome<-m1$model[,1]
    form <- as.character(m1$formula)
    form.2<-as.formula(paste(form[2],form[1],form[3],"|",form[3]))
    m1<-glm(m1$formula ,family="poisson",data=m1$data )
    m2<- suppressMessages(MASS::glm.nb(m1$formula ,data=m1$data ))
    m3 <- suppressMessages(pscl::zeroinfl(form.2,data=m1$data))
    m4 <- pscl::zeroinfl(form.2,data=m1$data, dist="negbin")
    bics<-c(BIC(m1),BIC(m2),BIC(m3),BIC(m4))
    aics<-c(AIC(m1),AIC(m2),AIC(m3),AIC(m4))
    out<-matrix(0,nrow=2*(ncol(m1$model)),ncol=16)
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
    zpoi<-round(apply(pscl::predprob(m3),2,mean)[1:length(y.range)],rounded)
    znb<-round(apply(pscl::predprob(m4),2,mean)[1:length(y.range)],rounded)
    poi<-obs-poi
    nb<-obs-nb
    zpoi<-obs-zpoi
    znb<-obs-znb
    type<-c(rep("Pois",length(y.range)),rep("NB",length(y.range)),rep("ZIP",length(y.range)),rep("ZNB",length(y.range)))
    dat<-data.frame(a=c(poi,nb,zpoi,znb),b=rep(y.range,4),type=type)

    if(use.color=="yes"){
      suppressMessages(suppressWarnings(
        out.plot<- ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x=b, y=a),group=type,linetype=type) +
          ggplot2::geom_hline(yintercept=0,color="gray80") +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position="bottom") +
          ggplot2::labs(x="Outcome",y="Observed - Predicted",color="") +
          ggplot2::scale_x_continuous(breaks=seq(min(y.range),max(y.range),round(length(y.range)/4)))  +
          ggplot2::geom_line(ggplot2::aes(linetype=type,color=type)) + ggplot2::labs(linetype="",color="") +
          ggplot2::scale_shape_manual(name="",values=c(22,21,23,24)) +
          ggplot2::theme(panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))))}else{
            suppressMessages(suppressWarnings(out.plot<- ggplot2::ggplot(data = dat, mapping = ggplot2::aes(x=b, y=a),group=type,linetype=type) +
                                                ggplot2::geom_hline(yintercept=0,color="gray80") +
                                                ggplot2::theme_bw() +
                                                ggplot2::theme(legend.position="bottom") +
                                                ggplot2::xlab("Outcome") +
                                                ggplot2::ylab("Observed - Predicted") +
                                                ggplot2::scale_x_continuous(breaks=seq(min(y.range),max(y.range),round(length(y.range)/4)))  +
                                                ggplot2::geom_line(ggplot2::aes(linetype=type)) + ggplot2::labs(linetype="") +
                                                ggplot2::scale_shape_manual(name="",values=c(22,21,23,24)) +
                                                ggplot2::theme(panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))))

          }

    models <- out[,c(1,2,5,6,9,10,13,14)]
    pdat <- data.frame(var=c("Intercept",names(m1$coef)[2:length(names(m1$coefficients))]),
                       coef=as.numeric(c(models[,1],models[,3],models[,5],models[,7])),
                       se=as.numeric(c(models[,2],models[,4],models[,6],models[,8])),
                       model=rep(c("Poisson","Neg Binom","ZIP","ZNB"),each=nrow(models)))

    pdat <- dplyr::mutate(pdat,max=coef+1.96*se,min=coef-1.96*se)

    pdat$var <- paste(rep(rep(c("c","z"),each=length(m1$coefficients)),4),pdat$var)
    pdat1<-pdat[1:(nrow(pdat)/4),]
    pdat1$var <- factor(pdat1$var,levels=paste(pdat1$var))
    pdat1[((nrow(pdat1)/2)+1):nrow(pdat1),2:3]<-NA
    p1<- suppressMessages(suppressWarnings(ggplot2::ggplot(data = pdat1,mapping = ggplot2::aes(y=var,x=coef,xmin=min,xmax=max)) +
                                             ggplot2::geom_vline(xintercept=0,color="gray80") +
                                             ggplot2::geom_pointrange() +
                                             ggplot2::theme_classic() + ggplot2::labs(x="",y="",title="Poisson")  + ggplot2::scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))))
    pdat2<-pdat[((nrow(pdat)/4)+1):((nrow(pdat)/2)),]
    pdat2$var <- factor(pdat2$var,levels=paste(pdat2$var))
    pdat2[((nrow(pdat2)/2)+1):nrow(pdat2),2:3]<-NA
    p2<- suppressMessages(suppressWarnings(ggplot2::ggplot(data = pdat2,mapping = ggplot2::aes(y=var,x=coef,xmin=min,xmax=max)) +
                                             ggplot2::geom_vline(xintercept=0,color="gray80") +
                                             ggplot2::geom_pointrange() +
                                             ggplot2::theme_classic() + ggplot2::labs(x="",y="",title="Neg Binom")  + ggplot2::scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))))
    pdat3<-pdat[((nrow(pdat)/2)+1):(((nrow(pdat)/4))*3),]
    pdat3$var <- factor(pdat3$var,levels=paste(pdat3$var))
    p3<- suppressMessages(suppressWarnings(ggplot2::ggplot(data = pdat3,mapping = ggplot2::aes(y=var,x=coef,xmin=min,xmax=max)) +
                                             ggplot2::geom_vline(xintercept=0,color="gray80") +
                                             ggplot2::geom_pointrange() +
                                             ggplot2::theme_classic() + ggplot2::labs(x="",y="",title="ZIP")  + ggplot2::scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))))
    pdat4<-pdat[(1+(((nrow(pdat)/4))*3)):nrow(pdat),]
    pdat4$var <- factor(pdat4$var,levels=paste(pdat4$var))
    p4<- suppressMessages(suppressWarnings(ggplot2::ggplot(data = pdat4,mapping = ggplot2::aes(y=var,x=coef,xmin=min,xmax=max)) +
                                             ggplot2::geom_vline(xintercept=0,color="gray80") +
                                             ggplot2::geom_pointrange() +
                                             ggplot2::theme_classic() + ggplot2::labs(x="",y="",title="ZINB")  + ggplot2::scale_x_continuous(limits=c(min(pdat$min),max(pdat$max)))))
    coef.pic<- suppressMessages(suppressWarnings(ggpubr::ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)))

    information.criterion<-rbind(bics,aics)
    rownames(information.criterion)<-c("BIC","AIC")
    colnames(information.criterion)<-c("Poisson","Neg Binom","ZIP","ZNB")

    outp<-list(ic=information.criterion,pic=out.plot,models=out,models.pic=coef.pic)} else{
      outp <- message("Model type should be Poisson regression as estimated by the glm function.")
    }
  return(outp)}















first.diff.fitted <- function (mod, design.matrix, compare, alpha = 0.05, rounded = 3,
                               bootstrap = "no", num.sample = 1000, prop.sample = 0.9, data,
                               seed = 1234,cum.probs="no")  {
  m1 <- NULL; design <- NULL
  if(is(mod,"lm") | is(mod,"glm") | is(mod,"multinom") | is(mod,"vglm")  |
     is(mod,"negbin") | is(mod,"zeroinfl") | is(mod,"hurdle") |
     is(mod,"glmerMod")  | is(mod,"clmm") | is(mod,"polr")){

    if (bootstrap == "no" & is(mod,"lm") | bootstrap == "no" & is(mod,"glm") |
        bootstrap == "no" & is(mod,"negbin") | bootstrap == "no" & is(mod,"lme")){

      des1 <- design.matrix[compare[1:2], ]
      f2 <- margins.dat(mod, design.matrix[compare[1:2], ])
      if (is(mod,"lm")) {
        f1 <- function(x) predict(x, type = "response", newdata = des1[1,
        ]) - predict(x, type = "response", newdata = des1[2,
        ])
      }
      else {
        f1 <- function(x) predict(x, type = "response", newdata = des1[1,
        ]) - predict(x, type = "response", newdata = des1[2,
        ])
      }
      out <- marginaleffects::hypotheses(mod,f1)
      out <- c(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[5] <- out[1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      out[6] <- out[1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      names(out) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- as.matrix(t(out))
      out <- round(out, rounded)
      out <- as.data.frame(out)
      if (sum(class(mod) == "nnet") == 1) {
        out <- cbind(out, levels(f2[nrow(design.matrix)/2,
                                    ncol(design.matrix) + 1]))
        colnames(out)[ncol(out)] <- "dv"
      }
      if (length(compare) > 2) {
        lc <- length(compare)/2
        coms <- seq(1, length(compare), 2)
        coms <- coms[-1]
        for (j in 1:(lc - 1)) {
          compare2 <- compare[coms[j]:(coms[j] + 1)]
          f2 <- margins.dat(mod, design.matrix[compare2,
          ])
          if (sum(class(mod) == "lm") > 0) {
            f1 <- function(x) predict(x, type = "response",
                                      newdata = design.matrix[compare2[1], ]) -
              predict(x, type = "response", newdata = design.matrix[compare2[2],
              ])
          }
          else {
            f1 <- function(x) predict(x, type = "response",
                                      newdata = design.matrix[compare2[1], ]) -
              predict(x, type = "response", newdata = design.matrix[compare2[2],
              ])
          }
          out2 <- marginaleffects::hypotheses(mod,f1)
          out2 <- c(out2$estimate,out2$std.error,out2$statistic,out2$p.value,out2$conf.low,out2$conf.high)
          out2[5] <- out2[1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[2]
          out2[6] <- out2[1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[2]
          names(out2) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
          out2 <- as.matrix(t(out2))
          out2 <- round(out2, rounded)
          out2 <- as.data.frame(out2)
          if (sum(class(mod) == "nnet") == 1) {
            out2 <- cbind(out2, levels(f2[nrow(design.matrix)/2,
                                          ncol(design.matrix) + 1]))
            colnames(out2)[ncol(out2)] <- "dv"
          }
          out <- rbind(out, out2)
        }
      }
    }else if (bootstrap == "no" & is(mod,"multinom")){
      des1 <- design.matrix[compare[1:2], ]
      f1 <- function(x) predict(x, type = "probs", newdata = des1[1,]) - predict(x, type = "probs", newdata = des1[2,])
      out <- marginaleffects::hypotheses(mod,f1)
      out <- data.frame(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[,5] <- out[,1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      out[,6] <- out[,1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      names(out) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- round(out, rounded)
      if (length(compare) > 2) {
        lc <- length(compare)/2
        coms <- seq(1, length(compare), 2)
        coms <- coms[-1]
        for (j in 1:(lc - 1)) {
          compare2 <- compare[coms[j]:(coms[j] + 1)]
          f1 <- function(x) predict(x, type = "probs",
                                    newdata = design.matrix[compare2[1], ]) -
            predict(x, type = "probs", newdata = design.matrix[compare2[2],
            ])
          out2 <- marginaleffects::hypotheses(mod,f1)
          out2 <- data.frame(out2$estimate,out2$std.error,out2$statistic,out2$p.value,out2$conf.low,out2$conf.high)
          out2[,5] <- out2[,1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[,2]
          out2[,6] <- out2[,1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[,2]
          names(out2) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
          out2 <- round(out2, rounded)
          out <- rbind(out, out2)
        }
      }
      out <- data.frame(dv=names(predict(mod,type = "probs", newdata = des1[1,])),out)
      names(out) <- c("Outcome Level","First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")

    } else if (bootstrap == "no" & is(mod,"vglm") | bootstrap == "no" & is(mod,"zeroinfl") | bootstrap == "no" & is(mod,"hurdle")){
      out<-("Model not suppported with Delta method. Use bootstrapping.")
    } else if (bootstrap == "no" & is(mod,"glmerMod")){
      f1 <- function(x) predict(x, type = "response", newdata = design.matrix[compare[1],],re.form=NA) - predict(x, type = "response", newdata = design.matrix[compare[2],],re.form=NA)
      out <- marginaleffects::hypotheses(mod,f1)
      out <- c(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[5] <- out[1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      out[6] <- out[1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      names(out) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- as.matrix(t(out))
      out <- round(out, rounded)
      out <- as.data.frame(out)
      if (length(compare) > 2) {
        lc <- length(compare)/2
        coms <- seq(1, length(compare), 2)
        coms <- coms[-1]
        for (j in 1:(lc - 1)) {
          compare2 <- compare[coms[j]:(coms[j] + 1)]
          f1 <- function(x) predict(x, type = "response", newdata = design.matrix[compare2[1],],re.form=NA) - predict(x, type = "response", newdata = design.matrix[compare2[2],],re.form=NA)

          out2 <- marginaleffects::hypotheses(mod,f1)
          out2 <- c(out2$estimate,out2$std.error,out2$statistic,out2$p.value,out2$conf.low,out2$conf.high)
          out2[5] <- out2[1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[2]
          out2[6] <- out2[1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[2]
          names(out2) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
          out2 <- as.matrix(t(out2))
          out2 <- round(out2, rounded)
          out2 <- as.data.frame(out2)
          out <- rbind(out, out2)}
      }
    }else if (bootstrap == "no" & is(mod,"polr") & cum.probs=="no"){
      des1 <- design.matrix[compare[1:2], ]
      f1 <- function(x) predict(x, type = "probs", newdata = des1[1,
      ]) - predict(x, type = "probs", newdata = des1[2,
      ])
      out <- marginaleffects::hypotheses(mod,f1)
      out <- data.frame(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[,5] <- out[,1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      out[,6] <- out[,1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      names(out) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- round(out, rounded)
      if (length(compare) > 2) {
        lc <- length(compare)/2
        coms <- seq(1, length(compare), 2)
        coms <- coms[-1]
        for (j in 1:(lc - 1)) {
          compare2 <- compare[coms[j]:(coms[j] + 1)]
          f1 <- function(x) predict(x, type = "probs", newdata = design.matrix[compare2[1],]) - predict(x, type = "probs", newdata = design.matrix[compare2[2],])
          out2 <- marginaleffects::hypotheses(mod,f1)
          out2 <- data.frame(out2$estimate,out2$std.error,out2$statistic,out2$p.value,out2$conf.low,out2$conf.high)
          out2[,5] <- out2[,1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[,2]
          out2[,6] <- out2[,1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
            out2[,2]
          names(out2) <- c("First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
          out2 <- round(out2, rounded)
          out <- rbind(out, out2)
        }
      }
      out <- data.frame(outcome.level=mod$lev,out)
      names(out) <- c("Outcome Level","First Difference","Standard Error","Statistic","p-value" ,"ll", "ul")

    }else if (bootstrap == "no" & is(mod,"polr") & cum.probs=="yes"){
      out <- message("Cumulative probabilities are not supported with parametric inference/delta method. Try bootstrapping or non-cumulative probabilities.")
    }
    else if (bootstrap == "no" & is(mod,"clmm")){out <- "Delta Method is not supported for mixed effects ordinal logistic regression. Use bootstrapping :("}
    else if (bootstrap == "yes" & is(mod,"glmerMod")){out<-"Bootstrapping is not supported. Use the delta method." }
    else if (bootstrap=="yes"){
      # Start Bootstrap code
      # Linear Regression
      if(is(mod,"lm")){
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "response", newdata = des1[1,]) -
          predict(mod, type = "response", newdata = des1[2,])
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- lm(formula(mod),data=fd.model2)
          obs.diffi<- predict(fd.modi, type = "response", newdata = des1[1,]) -
            predict(fd.modi, type = "response", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        fd.dist[,1] <- sort(fd.dist[,1])
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(sd(fd.dist),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"glm")){
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "response", newdata = des1[1,]) -
          predict(mod, type = "response", newdata = des1[2,])
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- glm(formula(mod),data=fd.model2,family=family(mod))
          obs.diffi<- predict(fd.modi, type = "response", newdata = des1[1,]) -
            predict(fd.modi, type = "response", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        fd.dist[,1] <- sort(fd.dist[,1])
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(sd(fd.dist),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"polr")){
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "response", newdata = des1[1,]) -
          predict(mod, type = "response", newdata = des1[2,])
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- MASS::polr(formula(mod),data=fd.model2, Hess=TRUE)
          obs.diffi<- predict(fd.modi, type = "response", newdata = des1[1,]) -
            predict(fd.modi, type = "response", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"multinom")){
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "probs", newdata = des1[1,]) -
          predict(mod, type = "probs", newdata = des1[2,])
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
          fd.modi <- nnet::multinom(formula(mod),data=fd.model2)
          obs.diffi<- predict(fd.modi, type = "probs", newdata = des1[1,]) -
            predict(fd.modi, type = "probs", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
        names(out) <- c("First Difference","SD of the distribution","ll.boot","ul.boot")
      } else if(is(mod,"vglm")){
        out <- ("Partial proportional odds model not suppported in a general way. An example is provided in the following script: [first.diff.fitted.partial.prop.odds.e.g.r]")
      } else if(is(mod,"negbin")){
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "response", newdata = des1[1,]) -
          predict(mod, type = "response", newdata = des1[2,])
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- MASS::glm.nb(formula(mod),data=fd.model2)
          obs.diffi<- predict(fd.modi, type = "response", newdata = des1[1,]) -
            predict(fd.modi, type = "response", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"zeroinfl")) {
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "response", newdata = des1[1,]) -
          predict(mod, type = "response", newdata = des1[2,])
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- pscl::zeroinfl(formula(mod),dist=mod$dist,data=fd.model2)
          obs.diffi<- predict(fd.modi, type = "response", newdata = des1[1,]) -
            predict(fd.modi, type = "response", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
        names(out) <- c("First Difference","SD of the distribution","ll.boot","ul.boot")
      }  else if(is(mod,"hurdle")) {
        des1 <- design.matrix[compare[1:2], ]
        obs.diff<- predict(mod, type = "response", newdata = des1[1,]) -
          predict(mod, type = "response", newdata = des1[2,])
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- pscl::hurdle(formula(mod),dist=mod$dist$count,data=fd.model2)
          obs.diffi<- predict(fd.modi, type = "response", newdata = des1[1,]) -
            predict(fd.modi, type = "response", newdata = des1[2,])
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
        names(out) <- c("First Difference","SD of the distribution","ll.boot","ul.boot")
      }  else if(is(mod,"clmm") & cum.probs=="no"){
        obs.diff<- suppressMessages(data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(design[compare[1],])))$prob -
          data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(design[compare[2],])))$prob)
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- ordinal::clmm(formula(mod), data=fd.model2,na.action=na.omit)
          obs.diffi<- suppressMessages(data.frame(emmeans::emmeans(fd.modi,~dv,mode="prob",at=as.list(design[compare[1],])))$prob -
            data.frame(emmeans::emmeans(fd.modi,~dv,mode="prob",at=as.list(design[compare[2],])))$prob)
          fd.dist[i,] <- obs.diffi}
        fd.dist[,1] <- sort(fd.dist[,1])
        out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(sd(fd.dist),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"clmm") & cum.probs=="yes"){
        obs.diff<- suppressMessages(data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(design[compare[1],])))$cumprob -
          data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(design[compare[2],])))$cumprob)
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- ordinal::clmm(formula(mod), data=fd.model2,na.action=na.omit)
          obs.diffi<- suppressMessages(data.frame(emmeans::emmeans(fd.modi,~cut,mode="cum.prob",at=as.list(design[compare[1],])))$cumprob -
            data.frame(emmeans::emmeans(fd.modi,~cut,mode="cum.prob",at=as.list(design[compare[2],])))$cumprob)
          fd.dist[i,] <- obs.diffi}
        fd.dist <- apply(fd.dist,2,FUN="sort")
        out <- data.frame(first.diff=round(obs.diff,rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      }
    }
    return(out)
  }else{message("Model type is not supported.")}}

second.diff.fitted <- function (mod, design.matrix, compare, alpha = 0.05, rounded = 3,
                                bootstrap = "no", num.sample = 1000, prop.sample = 0.9, data,
                                seed = 1234,cum.probs="no") {
  if(is(mod,"lm") | is(mod,"glm") | is(mod,"polr") | is(mod,"multinom") | is(mod,"vglm")  |
     is(mod,"negbin") | is(mod,"zeroinfl") |
     is(mod,"hurdle") | is(mod,"glmerMod")  | is(mod,"clmm")){

    if (bootstrap == "no" & is(mod,"lm") | bootstrap == "no" & is(mod,"glm") |
        bootstrap == "no" & is(mod,"negbin")){
      if (is(mod,"lm")) {
        f1 <- function(x) (predict(x, type = "response",
                                   newdata = design.matrix[compare[1], ]) - predict(x,
                                                                                    type = "response", newdata = design.matrix[compare[2],
                                                                                    ])) - (predict(x, type = "response", newdata = design.matrix[compare[3],
                                                                                    ]) - predict(x, type = "response", newdata = design.matrix[compare[4],
                                                                                    ]))
      }
      else {
        f1 <- function(x) (predict(x, type = "response", newdata = design.matrix[compare[1],
        ]) - predict(x, type = "response", newdata = design.matrix[compare[2],
        ])) - (predict(x, type = "response", newdata = design.matrix[compare[3],
        ]) - predict(x, type = "response", newdata = design.matrix[compare[4],
        ]))
      }
      out <- marginaleffects::hypotheses(mod,f1)
      out <- c(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[5] <- out[1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      out[6] <- out[1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      names(out) <- c("Second Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- as.matrix(t(out))
      out <- round(out, rounded)
      out <- as.data.frame(out)
    } else if (bootstrap == "no" & is(mod,"multinom")){
      f1 <- function(x) (predict(x, type = "probs", newdata = design.matrix[compare[1],]) -
                           predict(x, type = "probs", newdata = design.matrix[compare[2],])) -
        (predict(x, type = "probs", newdata = design.matrix[compare[3],]) -
           predict(x, type = "probs", newdata = design.matrix[compare[4],]))
      out <- marginaleffects::hypotheses(mod,f1)
      out <- data.frame(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[,5] <- out[,1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      out[,6] <- out[,1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      names(out) <- c("Second Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- round(out, rounded)
      out <- data.frame(dv=names(predict(mod,type = "probs", newdata = design.matrix[1,])),out)
      names(out) <- c("Outcome Level","Second Difference","Standard Error","Statistic","p-value" ,"ll", "ul")

    } else if (bootstrap == "no" & is(mod,"glmerMod")){
      f1 <- function(x) (predict(mod, type = "response",newdata = design.matrix[compare[1], ],re.form=NA) -
                           predict(mod,type = "response", newdata = design.matrix[compare[2],],re.form=NA)) -
        (predict(mod, type = "response", newdata = design.matrix[compare[3],],re.form=NA) -
           predict(x, type = "response", newdata = design.matrix[compare[4], ],re.form=NA))
      out <- marginaleffects::hypotheses(mod,f1)
      out <- c(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[5] <- out[1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      out[6] <- out[1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[2]
      names(out) <- c("Second Difference","Standard Error","Statistic","p-value" ,"ll", "ul")
      out <- as.matrix(t(out))
      out <- round(out, rounded)
      out <- as.data.frame(out)
    }else if (bootstrap == "no" & is(mod,"polr") & cum.probs=="no"){
      des1 <- design.matrix[compare[1:4], ]
      f1 <- function(x) (predict(x, type = "probs", newdata = des1[1,]) - predict(x, type = "probs", newdata = des1[2,])) -
        (predict(x, type = "probs", newdata = des1[3,]) - predict(x, type = "probs", newdata = des1[4,]))
      out <- marginaleffects::hypotheses(mod,f1)
      out <- data.frame(out$estimate,out$std.error,out$statistic,out$p.value,out$conf.low,out$conf.high)
      out[,5] <- out[,1] - qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      out[,6] <- out[,1] + qnorm(1 - (alpha/2), lower.tail = TRUE) *
        out[,2]
      out <- round(out, rounded)
      out <- data.frame(outcome.level=mod$lev,out)
      names(out) <- c("Outcome Level","Second Difference","Standard Error","Statistic","p-value" ,"ll", "ul")

    }else if (bootstrap == "no" & is(mod,"polr") & cum.probs=="yes"){
      out <- message("Cumulative probabilities are not supported with parametric inference/delta method. Try bootstrapping or non-cumulative probabilities.")
    }else if (bootstrap == "no" & is(mod,"vglm") | bootstrap == "no" & is(mod,"zeroinfl") | bootstrap == "no" & is(mod,"hurdle") | bootstrap == "no" & is(mod,"clmm")){
      out<-message("Model not suppported with Delta method. Use bootstrapping.")
    } else if (bootstrap=="yes"){

      # Start Bootstrap code
      # Linear Regression
      if(is(mod,"lm")){

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "response", newdata = des1[1,]) -
                      predict(mod, type = "response", newdata = des1[2,])) - (predict(mod, type = "response", newdata = des1[3,]) -
                                                                                predict(mod, type = "response", newdata = des1[4,]))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- lm(formula(mod),data=fd.model2)
          obs.diffi<- (predict(fd.modi, type = "response", newdata = des1[1,]) -
                         predict(fd.modi, type = "response", newdata = des1[2,])) - (predict(fd.modi, type = "response", newdata = des1[3,]) -
                                                                                       predict(fd.modi, type = "response", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        fd.dist[,1] <- sort(fd.dist[,1])
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(sd(fd.dist),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"glm")){

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "response", newdata = des1[1,]) -
                      predict(mod, type = "response", newdata = des1[2,])) - (predict(mod, type = "response", newdata = des1[3,]) -
                                                                                predict(mod, type = "response", newdata = des1[4,]))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- glm(formula(mod),data=fd.model2,family=family(mod))
          obs.diffi<- (predict(fd.modi, type = "response", newdata = des1[1,]) -
                         predict(fd.modi, type = "response", newdata = des1[2,])) - (predict(fd.modi, type = "response", newdata = des1[3,]) -
                                                                                       predict(fd.modi, type = "response", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        fd.dist[,1] <- sort(fd.dist[,1])
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(sd(fd.dist),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"polr")){

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "response", newdata = des1[1,]) -
                      predict(mod, type = "response", newdata = des1[2,])) - (predict(mod, type = "response", newdata = des1[3,]) -
                                                                                predict(mod, type = "response", newdata = des1[4,]))

        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- MASS::polr(formula(mod),data=fd.model2, Hess=TRUE)
          obs.diffi<- (predict(fd.modi, type = "response", newdata = des1[1,]) -
                         predict(fd.modi, type = "response", newdata = des1[2,])) - (predict(fd.modi, type = "response", newdata = des1[3,]) -
                                                                                       predict(fd.modi, type = "response", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"multinom")){

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "probs", newdata = des1[1,]) -
                      predict(mod, type = "probs", newdata = des1[2,])) - (predict(mod, type = "probs", newdata = des1[3,]) -
                                                                             predict(mod, type = "probs", newdata = des1[4,]))
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- data[sample(1:nrow(data),round(prop.sample*nrow(data),0),replace=TRUE),]
          fd.modi <- nnet::multinom(formula(mod),data=fd.model2)
          obs.diffi<- (predict(fd.modi, type = "probs", newdata = des1[1,]) -
                         predict(fd.modi, type = "probs", newdata = des1[2,])) - (predict(fd.modi, type = "probs", newdata = des1[3,]) -
                                                                                    predict(fd.modi, type = "probs", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
        names(out) <- c("Second Difference","SD of the distribution","ll.boot","ul.boot")
      } else if(is(mod,"vglm")){

        out <- ("Partial proportional odds model not suppported in a general way. An example is provided in the following script: [first.diff.fitted.partial.prop.odds.e.g.r]")
      } else if(is(mod,"negbin")){

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "response", newdata = des1[1,]) -
                      predict(mod, type = "response", newdata = des1[2,])) - (predict(mod, type = "response", newdata = des1[3,]) -
                                                                                predict(mod, type = "response", newdata = des1[4,]))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- MASS::glm.nb(formula(mod),data=fd.model2)
          obs.diffi<- (predict(fd.modi, type = "response", newdata = des1[1,]) -
                         predict(fd.modi, type = "response", newdata = des1[2,])) - (predict(fd.modi, type = "response", newdata = des1[3,]) -
                                                                                       predict(fd.modi, type = "response", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"zeroinfl")) {

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "response", newdata = des1[1,]) -
                      predict(mod, type = "response", newdata = des1[2,])) - (predict(mod, type = "response", newdata = des1[3,]) -
                                                                                predict(mod, type = "response", newdata = des1[4,]))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <-  pscl::zeroinfl(formula(mod),dist=mod$dist,data=fd.model2)
          obs.diffi<- (predict(fd.modi, type = "response", newdata = des1[1,]) -
                         predict(fd.modi, type = "response", newdata = des1[2,])) - (predict(fd.modi, type = "response", newdata = des1[3,]) -
                                                                                       predict(fd.modi, type = "response", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
        names(out) <- c("Second Difference","SD of the distribution","ll.boot","ul.boot")
      } else if(is(mod,"hurdle")) {

        des1 <- design.matrix[compare[1:4], ]
        obs.diff<- (predict(mod, type = "response", newdata = des1[1,]) -
                      predict(mod, type = "response", newdata = des1[2,])) - (predict(mod, type = "response", newdata = des1[3,]) -
                                                                                predict(mod, type = "response", newdata = des1[4,]))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- pscl::hurdle(formula(mod),dist=mod$dist$count,data=fd.model2)
          obs.diffi<- (predict(fd.modi, type = "response", newdata = des1[1,]) -
                         predict(fd.modi, type = "response", newdata = des1[2,])) - (predict(fd.modi, type = "response", newdata = des1[3,]) -
                                                                                       predict(fd.modi, type = "response", newdata = des1[4,]))
          fd.dist[i,] <- obs.diffi}
        for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])}
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
        names(out) <- c("Second Difference","SD of the distribution","ll.boot","ul.boot")
      }else if(is(mod,"clmm") & cum.probs=="no"){
        obs.diff<- suppressMessages((data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(design.matrix[compare[1],])))$prob -
                      data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(design.matrix[compare[2],])))$prob) -
          (data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(design.matrix[compare[3],])))$prob -
             data.frame(emmeans::emmeans(mod,~dv,mode="prob",at=as.list(design.matrix[compare[4],])))$prob))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- ordinal::clmm(formula(mod), data=fd.model2,na.action=na.omit)
          obs.diffi<- suppressMessages((data.frame(emmeans::emmeans(fd.modi,~dv,mode="prob",at=as.list(design.matrix[compare[1],])))$prob -
                         data.frame(emmeans::emmeans(fd.modi,~dv,mode="prob",at=as.list(design.matrix[compare[2],])))$prob) -
            (data.frame(emmeans::emmeans(fd.modi,~dv,mode="prob",at=as.list(design.matrix[compare[3],])))$prob -
               data.frame(emmeans::emmeans(fd.modi,~dv,mode="prob",at=as.list(design.matrix[compare[4],])))$prob))
          fd.dist[i,] <- obs.diffi}
        fd.dist <- apply(fd.dist,2,FUN="sort")
        out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(sd(fd.dist),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      } else if(is(mod,"clmm") & cum.probs=="yes"){
        obs.diff<- suppressMessages((data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(design.matrix[compare[1],])))$cumprob -
                      data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(design.matrix[compare[2],])))$cumprob) -
          (data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(design.matrix[compare[3],])))$cumprob -
             data.frame(emmeans::emmeans(mod,~cut,mode="cum.prob",at=as.list(design.matrix[compare[4],])))$cumprob))
        fd.model<-mod$model
        fd.dist <-matrix(NA,nrow=num.sample,ncol=length(obs.diff))
        for(i in 1:num.sample){
          set.seed(seed + i);  fd.model2 <- fd.model[sample(1:nrow(fd.model),round(prop.sample*nrow(fd.model),0),replace=TRUE),]
          fd.modi <- ordinal::clmm(formula(mod), data=fd.model2,na.action=na.omit)
          obs.diffi<- suppressMessages((data.frame(emmeans::emmeans(fd.modi,~cut,mode="cum.prob",at=as.list(design.matrix[compare[1],])))$cumprob -
                         data.frame(emmeans::emmeans(fd.modi,~cut,mode="cum.prob",at=as.list(design.matrix[compare[2],])))$cumprob) -
            (data.frame(emmeans::emmeans(fd.modi,~cut,mode="cum.prob",at=as.list(design.matrix[compare[3],])))$cumprob -
               data.frame(emmeans::emmeans(fd.modi,~cut,mode="cum.prob",at=as.list(design.matrix[compare[4],])))$cumprob))
          fd.dist[i,] <- obs.diffi}
        fd.dist <- apply(fd.dist,2,FUN="sort")
        out <- data.frame(second.diff=round(obs.diff,rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha/2)),],rounded))
      }else if(is(mod,"glmerMod")) {
        out <- "Model type not supported with bootstrapping. Use Delta method."
      }
    }
    return(out)
  }else{message("Model type is not supported.")}
}

compare.margins <- function(margins,margins.ses,seed=1234,rounded=3,nsim=10000){
  difference <- margins[1] - margins[2]
  if(difference>0){
    set.seed(seed); p.value<-sum((rnorm(nsim,mean=margins[1],sd=margins.ses[1]) - rnorm(nsim,mean=margins[2],sd=margins.ses[2])) < 0) /nsim
  }else{
    set.seed(seed); p.value<-sum((rnorm(nsim,mean=margins[1],sd=margins.ses[1]) - rnorm(nsim,mean=margins[2],sd=margins.ses[2])) > 0) /nsim
  }
  out <- c(difference,p.value)
  names(out) <- c("Difference","p-value")
  out <- as.data.frame(t(out))
  return(out)
}

rubins.rule <- function(std.errors){
  r.r.std.error<-sqrt(mean(std.errors^2) + var(std.errors) + var(std.errors)/length(std.errors))
  return(r.r.std.error)}

diagn <- function(model){
  if(is(model,"glm")){
    pearsonres<-residuals.glm(model,type="pearson")
    h<-data.frame(influence.measures(model)$infmat)$hat
    stdpres<-pearsonres/(sqrt(1-h))
    deltabeta<-(pearsonres^2*h)/(1-h)^2
    deltabeta <- deltabeta/length(coef(model))
    obs = 1:length(pearsonres)
    devres<-residuals.glm(model,type="deviance")
    out=data.frame(pearsonres,h,stdpres,deltabeta,obs,devres)
  }
  else if(is(model,"negbin")){
    pearsonres<-residuals.glm(model,type="pearson")
    h<-data.frame(influence.measures(model)$infmat)$hat
    stdpres<-pearsonres/(sqrt(1-h))
    deltabeta<-(pearsonres^2*h)/(1-h)^2
    deltabeta <- deltabeta/length(coef(model))
    obs = 1:length(pearsonres)
    devres<-residuals.glm(model,type="deviance")
    out=data.frame(pearsonres,h,stdpres,deltabeta,obs,devres)
  }

  else if(is(model,"polr")){
    out=c("Ordinal model: you should binarize it for diagnostics.")
  }
  else if(is(model,"multinom")){
    out=c("Multinomial model: you should binarize it for diagnostics.")
  }
  else if(is(model,"vglm")){
    out=c("Partial Proportional Odds model: you should binarize it for diagnostics.")
  }
  else if(is(model,"zeroinfl")){
    pearsonres<-residuals(model,type="pearson")
    obs = 1:length(pearsonres)
    out=data.frame(pearsonres,obs)
  }
  else if(is(model,"hurdle")){
    pearsonres<-residuals(model,type="pearson")
    obs = 1:length(pearsonres)
    out=data.frame(pearsonres,obs)
  }
  else{out="Model type not supported."}
  return(out)}

margins.dat.clogit <-function (mod, design.matrix, run.boot = "no", num.sample = 1000,
                               prop.sample = 0.9, alpha = 0.05, seed = 1234, rounded = 3) {
  '(strata)' <- NULL
  if(is(mod,"clogistic")){
    coefs <- as.numeric(na.omit(coef(mod)))
    des <- dplyr::mutate(design.matrix, lp = exp(as.matrix(design.matrix) %*%
                                            coefs), probs = lp/sum(lp))
    sims <- matrix(NA, ncol = nrow(design.matrix), nrow = num.sample)
    vcovc <- vcov(mod)
    sl <- which(is.na(coef(mod)))
    if (length(sl) > 0) {
      vcovc <- vcovc[-sl, -sl]
    }
    for (i in 1:num.sample) {
      set.seed(seed + i)
      lp <- exp(as.matrix(design.matrix) %*% MASS::mvrnorm(mu = coefs,
                                                     Sigma = vcovc))
      sims[i, ] <- lp/sum(lp)
    }
    sims <- apply(sims, 2, FUN = "sort")
    des <- dplyr::mutate(des, ll = sims[alpha/2 * nrow(sims), ], ul = sims[(1 -
                                                                       alpha/2) * nrow(sims), ], se = apply(sims, 2, FUN = "sd"))
    out <- round(des, rounded)
    if (run.boot == "yes") {
      boot.dist <- matrix(NA, nrow = num.sample, ncol = nrow(des))
      for (i in 1:num.sample) {
        set.seed(seed + i)
        mod2 <- mod$model[sample(1:nrow(mod$model), round(prop.sample *
                                                            nrow(mod$model), 0), replace = FALSE), ]
        m.1 <- Epi::clogistic(mod$formula, strata = `(strata)`,
                         data = mod2)
        coefs2 <- as.numeric(na.omit(coef(m.1)))
        design20 <- dplyr::mutate(des, lp = exp(as.matrix(design.matrix) %*%
                                           coefs2), probs = lp/sum(lp))
        boot.dist[i, ] <- design20[, ncol(design20)]
      }
      boot.dist[, 1] <- sort(boot.dist[, 1])
      if (ncol(boot.dist) > 1) {
        for (i in 2:ncol(boot.dist)) {
          boot.dist[, i] <- sort(boot.dist[, i])
        }
      }
      lower.limit <- boot.dist[nrow(boot.dist) * (alpha/2),
      ]
      upper.limit <- boot.dist[nrow(boot.dist) * (1 - (alpha/2)),
      ]
      des <- round(dplyr::mutate(des, ll.boot = lower.limit, ul.boot = upper.limit),
                   rounded)
      out <- list(des = des, boot.dist = boot.dist)
    }
    return(out)} else if(is(mod,"mlogit")){
      des <- dplyr::mutate(design.matrix, probs= predict(mod,newdata=design.matrix))
      m2 <- mod
      sims <- matrix(NA, ncol = nrow(design.matrix), nrow = num.sample)
      vcovc <- vcov(mod)
      for (i in 1:num.sample) {
        set.seed(seed + i)
        m2$coefficients <- MASS::mvrnorm(mu = coef(mod), Sigma = vcov(mod))
        sims[i, ] <- predict(m2,newdata=design.matrix)}
      sims <- apply(sims, 2, FUN = "sort")
      des <- dplyr::mutate(des, ll = sims[alpha/2 * nrow(sims), ], ul = sims[(1 -
                                                                         alpha/2) * nrow(sims), ], se = apply(sims, 2, FUN = "sd"))
      out <- round(des, rounded)
      return(out)
    } else{message("Model type is not supported.")}}






