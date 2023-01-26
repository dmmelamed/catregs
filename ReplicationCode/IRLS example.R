
# Code for iteratively reweighted least squares estimation of logistic regression coefficients
# Begin with the Baseline logistic regression model from Chapter 5. Up to line 24 reproduces the model using the glm function
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


y <- dat$safe
X2 <- select(dat,religious,minority,female,age,emp1,emp2)
X2 <- data.frame(constant=1,X2)
X2 <- as.matrix(X2)
# IRLS starts with OLS estimates
b<-solve(t(X2)  %*% X2)%*%t(X2)%*%y


linear.predictors <- X2 %*% b # Compute linear predictors
phat<-1/(1+exp(-1*(as.vector(linear.predictors)))) # Take logistic transformation
V<-diag((phat*(1-phat))) # Generate weight matrix
z<-X2%*%b+solve(V)%*%(y-phat) # Define pseudovalues (see Hosmer and Lemeshow)
b2<-solve(t(X2) %*% V %*% X2)%*%t(X2)%*%V%*%z # Update the coefficients
cbind(b,b2) # Put the LPM coefficients next to the new ones
b<-b2 # Update the coefficients
# Lines 32-41 is an iteration of reweighted least squares (well, just "weighted" least squares at this point)

# Highlight and run lines 45-51, and repeat until the coefficients no longer change. It takes 5 iterations for the estimates to converge to all decimal places shown
linear.predictors <- X2 %*% b
phat<-1/(1+exp(-1*(as.vector(linear.predictors))))
V<-diag((phat*(1-phat)))
z<-X2%*%b+solve(V)%*%(y-phat)
b2<-solve(t(X2) %*% V %*% X2)%*%t(X2)%*%V%*%z
cbind(b,b2)
b<-b2

cbind(b,coef(m1)) # The IRLS estimates are identical to the estimates derived from glm

# Standard Errors
vcov<-sqrt(diag(solve(t(X2) %*% V %*% X2)))
# As a check, compare what we just computed to the software estimates of the coefficients
cbind(vcov,sqrt(diag(vcov(m1))))
