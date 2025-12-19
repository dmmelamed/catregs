# Example script for differences in fitted values from a mixed effects multinomial logistic regression model

rm(list=ls()) #clear the workspace
require(catregs) # load required packages
require(tidyverse)

data("ess")
ess <- mutate(ess,important.matters.people=factor(important.matters.people))
ess$important.matters.people[which(ess$important.matters.people=="Don't know")] <- NA
ess$important.matters.people[which(ess$important.matters.people=="Refusal")] <- NA
ess$important.matters.people <- factor(ess$important.matters.people,
                                         levels = c("None","1","2","3","4-6","7-9","10 or more"))
ess <- mutate(ess,dv=important.matters.people)
ess <- mutate(ess,male=ifelse(gender=="Male",1,0),
              min=ifelse(ethnic.minority=="Yes",1,0))
ess <- mutate(ess,mc = male*num.children)
ess <- select(ess,dv,male,age,num.children,mc,min,education,country)
ess <- na.omit(ess)
mod <- mclogit::mblogit(dv ~ male + age + num.children + min + education, random = ~ 1 | country,data=ess)
summary(mod)

# E.g., look at first differences by male
table(ess$male)
design <- margins.des(mod,expand.grid(male=c(0,1)),data=ess) # Generate the design matxix

g <- as.matrix(summary(emmeans::emmeans(mod,specs="dv", at = as.list(design[1,]), weights = "proportional", mode = "prob")))
h <- as.matrix(summary(emmeans::emmeans(mod,specs="dv", at = as.list(design[2,]), weights = "proportional", mode = "prob")))
obs.diff <- as.numeric(g[,2]) - as.numeric(h[,2])
names(obs.diff) <- g[,1]
obs.diff


num.sample <- 100 # Set the number of samples to draw for the bootstrap
fd.dist <-matrix(NA,nr=num.sample,nc=length(obs.diff)) # An object that will become the bootstrapped distribution
prop.sample <- .9
for(i in 1:num.sample){ # Loop the bootstrapped distribution
  set.seed(1234 + i);  fd.model2 <- ess[sample(1:nrow(ess),round(prop.sample*nrow(ess),0),replace=TRUE),] # Create a sample from the data. Each sample is 90% of the original, with replacement
  fd.modi <-mclogit::mblogit(dv ~ male + age + num.children + mc + min + education, random = ~ 1 | country,data=fd.model2) # Estimate the same model as above, but with the reduced sample
  gi <- as.matrix(summary(emmeans::emmeans(fd.modi,specs="dv", at = as.list(design[1,]), weights = "proportional", mode = "prob")))
  hi <- as.matrix(summary(emmeans::emmeans(fd.modi,specs="dv", at = as.list(design[2,]), weights = "proportional", mode = "prob")))
  obs.diffi <- as.numeric(gi[,2]) - as.numeric(hi[,2])
  names(obs.diffi) <- gi[,1]
  fd.dist[i,] <- obs.diffi} # Save the differences from this pass of the loop

for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])} # Sort the bootstrapped distributions

alpha.level <- .05 # Set your alpha value for confidence intervals
rounded <- 3 # Set the number of decimals to round

out <- data.frame(first.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha.level/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha.level/2)),],rounded))
out # the first differences, the standard deviation of the bootstrapped distribution, and the lower and upper confidence limits.

###
# Above is done, just need to add/adjust the comments. Below I need to adapt this for second differences and add comments...
###


# Now for second differences. Respecify mod, but include an interaction
mod <- mclogit::mblogit(dv ~ male + num.children + age +  mc + min + education, random = ~ 1 | country,data=ess)
design <- margins.des(mod,expand.grid(male=c(0,1),num.children=0:5),data=ess) # Generate the design matxix

k <- as.matrix(summary(emmeans::emmeans(mod,specs="dv", at = as.list(design[1,]), weights = "proportional", mode = "prob")))
l <- as.matrix(summary(emmeans::emmeans(mod,specs="dv", at = as.list(design[2,]), weights = "proportional", mode = "prob")))
m <- as.matrix(summary(emmeans::emmeans(mod,specs="dv", at = as.list(design[3,]), weights = "proportional", mode = "prob")))
n <- as.matrix(summary(emmeans::emmeans(mod,specs="dv", at = as.list(design[4,]), weights = "proportional", mode = "prob")))
obs.diff <- (as.numeric(k[,2]) - as.numeric(l[,2])) - (as.numeric(m[,2]) - as.numeric(n[,2]))
names(obs.diff) <- k[,1]
obs.diff

num.sample <- 100 # Set the number of samples to draw for the bootstrap
fd.dist <-matrix(NA,nr=num.sample,nc=length(obs.diff)) # An object that will become the bootstrapped distribution
prop.sample <- .9
for(i in 1:num.sample){ # Loop the bootstrapped distribution
  set.seed(1234 + i);  fd.model2 <- ess[sample(1:nrow(ess),round(prop.sample*nrow(ess),0),replace=TRUE),] # Create a sample from the data. Each sample is 90% of the original, with replacement
  fd.modi <-mclogit::mblogit(dv ~ male + age + num.children + mc + min + education, random = ~ 1 | country,data=fd.model2) # Estimate the same model as above, but with the reduced sample
  ki <- as.matrix(summary(emmeans::emmeans(fd.modi,specs="dv", at = as.list(design[1,]), weights = "proportional", mode = "prob")))
  li <- as.matrix(summary(emmeans::emmeans(fd.modi,specs="dv", at = as.list(design[2,]), weights = "proportional", mode = "prob")))
  mi <- as.matrix(summary(emmeans::emmeans(fd.modi,specs="dv", at = as.list(design[3,]), weights = "proportional", mode = "prob")))
  ni <- as.matrix(summary(emmeans::emmeans(fd.modi,specs="dv", at = as.list(design[4,]), weights = "proportional", mode = "prob")))
  obs.diffi <- (as.numeric(ki[,2]) - as.numeric(li[,2])) - (as.numeric(mi[,2]) - as.numeric(ni[,2]))
  names(obs.diffi) <- ki[,1]
  fd.dist[i,] <- obs.diffi} # Save the differences from this pass of the loop

for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])} # Sort the bootstrapped distributions

alpha.level <- .05 # Set your alpha value for confidence intervals
rounded <- 4 # Set the number of decimals to round

out <- data.frame(second.diff=round(obs.diff,rounded),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha.level/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha.level/2)),],rounded))
out # the second differences, the standard deviation of the bootstrapped distribution, and the lower and upper confidence limits.


