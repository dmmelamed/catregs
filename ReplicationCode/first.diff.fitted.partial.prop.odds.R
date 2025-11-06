# Example script for differences in fitted values from a Partial Proportional odds model

rm(list=ls()) #clear the workspace
require(catregs) # load required packages
require(tidyverse)

#setwd("~/Desktop/School/Research/InPrep/Categorical Regression Book/Replicate Long and Freese/lf2")
#X <- read.csv("oldwarm2.csv")
# Import the "oldwarm2.csv" file and call it "X"

X <- mutate(X,warm=as.factor(warm))
require(VGAM)
m6 <- vglm(warm ~ yr89 + male + white + age + ed + prst,
           cumulative(parallel = FALSE ~ yr89 + male + age, reverse = FALSE),
           data = X) # Specify a partial proportional odds model
summary(m6)
# E.g., look at first differences by male
table(X$male)
require(nnet) # Specify a multinomial model to generate the design matrix
m5<-multinom(warm~ yr89 + male + white + age + ed + prst, data=X)

design <- margins.des(m5,expand.grid(male=c(0,1)),data=X) # Generate the design matxix



obs.diff<- predict(m6, type = "response", newdata = design[1,]) - 
  predict(m6, type = "response", newdata = design[2,]) # Compute the observed first difference

num.sample <- 1000 # Set the number of samples to draw for the bootstrap
fd.dist <-matrix(NA,nr=num.sample,nc=length(obs.diff)) # An object that will become the bootstrapped distribution
for(i in 1:num.sample){ # Loop the bootstrapped distribution
  set.seed(1982 + i);  fd.model2 <- X[sample(1:nrow(X),round(.9*nrow(X),0),replace=TRUE),] # Create a sample from the data. Each sample is 90% of the original, with replacement
  fd.modi <-vglm(warm ~ yr89 + male + white + age + ed + prst,
                 cumulative(parallel = FALSE ~ yr89 + male + age, reverse = FALSE),
                 data = fd.model2) # Estimate the same model as above, but with the reduced sample
  obs.diffi<- predict(fd.modi, type = "response", newdata = design[1,]) - 
    predict(fd.modi, type = "response", newdata = design[2,]) # COmpute and retain the first differences for the model in this pass of the loop
  fd.dist[i,] <- obs.diffi} # Save the differences from this pass of the loop
for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])} # Sort the bootstrapped distributions

alpha.level <- .05 # Set your alpha value for confidence intervals
rounded <- 3 # Set the number of decimals to round

out <- data.frame(first.diff=t(round(obs.diff,rounded)),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha.level/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha.level/2)),],rounded))
colnames(out)[1] <- "first.difference"
out # the first differences, the standard deviation of the bootstrapped distribution, and the lower and upper confidence limits.


# Now for second differences. Respecify m6, but include an interaction
m6 <- vglm(warm ~ yr89 + male*white + age + ed + prst,
           cumulative(parallel = FALSE ~ yr89 + male + age, reverse = FALSE),
           data = X)

design <- margins.des(m5,expand.grid(male=c(0,1),white=c(0,1)),data=X) # Design matrix with at least 4 rows...
obs.diff<- (predict(m6, type = "response", newdata = design[1,]) - 
  predict(m6, type = "response", newdata = design[2,])) - (predict(m6, type = "response", newdata = design[3,]) - 
  predict(m6, type = "response", newdata = design[4,])) # Compute the observed second difference

num.sample <- 500 # Number of samples to draw for the bootstrap
fd.dist <-matrix(NA,nr=num.sample,nc=length(obs.diff))  # An object that will become the bootstrapped distribution
for(i in 1:num.sample){  # Loop the bootstrapped distribution
  set.seed(1982 + i);  fd.model2 <- X[sample(1:nrow(X),round(.9*nrow(X),0),replace=TRUE),] # Create a sample from the data. Each sample is 90% of the original, with replacement
  fd.modi <-vglm(warm ~ yr89 + male*white + age + ed + prst,
                 cumulative(parallel = FALSE ~ yr89 + male + age, reverse = FALSE),
                 data = fd.model2) # Estimate the same model as above, but with the reduced sample
  obs.diffi<- (predict(fd.modi, type = "response", newdata = design[1,]) - 
    predict(fd.modi, type = "response", newdata = design[2,])) - (predict(fd.modi, type = "response", newdata = design[3,]) - 
                                                                    predict(fd.modi, type = "response", newdata = design[4,])) # COmpute and retain the second differences for the model in this pass of the loop
  fd.dist[i,] <- obs.diffi}  # Save the differences from this pass of the loop
for(i in 1:ncol(fd.dist)){fd.dist[,i] <- sort(fd.dist[,i])} # Sort the bootstrapped distributions

alpha.level <- .05
rounded <- 3
out <- data.frame(first.diff=t(round(obs.diff,rounded)),sd.boot.dist=round(apply(fd.dist,2,"sd"),rounded),ll.boot=round(fd.dist[nrow(fd.dist)*(alpha.level/2),],rounded),ul.boot=round(fd.dist[nrow(fd.dist)*(1-(alpha.level/2)),],rounded))
colnames(out)[1] <- "second.difference"
out # the second differences, the standard deviation of the bootstrapped distribution, and the lower and upper confidence limits.








