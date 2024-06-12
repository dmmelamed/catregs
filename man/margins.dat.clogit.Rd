\name{margins.dat.clogit}
\alias{margins.dat.clogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes predicted probabilities for conditional and rank-order/exploded logistic regression models. Inference is based upon simulation techniques (requires the MASS package). Alternatively, bootstrapping is an option for conditional logistic regression models.
}
\description{
Given a model object and a design matrix, this creates a data.frame of the design matrix, with predicted probabilities for each response category. Inferential information about the predict probabilities is supported with simulation. Bootstrapping may be added as an option for conditional logistic regression models.
}
\usage{
margins.dat.clogit(mod,design.matrix,run.boot="no",num.sample=1000,
prop.sample=.9,alpha=.05,seed=1234,rounded=3)}
\arguments{
  \item{mod}{A conditional logistic regression model as estimated in the Epi package or an Exploded logistic regression model as estimated in the mlogit package.
}
  \item{design.matrix}{Design matrix of values for the independent variables in the regression model. Unlike the design matrices in the margins.des function, the design matrix for a conditional logistic regression entails multiple rows, corresponding to the number of response options.}
  \item{run.boot}{Whether to compute confidence intervals around the predicted probabilities using bootstrapping. Defaul is "no."}
  \item{num.sample}{num.sample is the number samples drawn to compute the sampling distibution.}
  \item{prop.sample}{prop.sample is the proportion of the original sample to include in the sampling distibution samples. Default is .9}
    \item{alpha}{The alpha value for confidence intervals. Default is .05.}
    \item{seed}{Sets a seed so that random results are reprodicible.}
    \item{rounded}{How many decimal places to show in the output.}

}
\value{
  \item{des}{Returns a data.frame containing the design matrix and additional columns for the predicted probabilities.}
  \item{boot.dist}{The full bootstrapped distribution for the probabilities.}
}
\author{
David Melamed
}
\examples{
data("LF06travel")
m1 <- Epi::clogistic(choice ~ train + bus +
time  + invc, strata=id, data=LF06travel)
design <- data.frame(train=c(0,0,1),bus=c(0,1,0),time=200,invc=20)
design
margins.dat.clogit(m1,design)
}
