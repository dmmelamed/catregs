\name{margins.dat.clogit}
\alias{margins.dat.clogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes predicted probabilities for conditional logistic regression models, with or without bootstrapped sampling distributions.
}
\description{
Given a model object and a design matrix, this creates a data.frame of the design matrix, with predicted probabilities for each response category. Inferential information about the predict probabilities is supported with bootstrapping.
}
\usage{
margins.dat.clogit(mod,design.matrix,run.boot="no",num.sample=1000,prop.sample=.9,alpha=.05)}
\arguments{
  \item{mod}{A conditional logistic regression model, as estimated in the Epi package.
}
  \item{design.matrix}{Design matrix of values for the independent variables in the regression model. Unlike the design matrices in the margins.des function, the design matrix for a conditional logistic regression entails multiple rows, corresponding to the number of response options.}
  \item{run.boot}{Whether to compute confidence intervals around the predicted probabilities using bootstrapping. Defaul is "no."}
  \item{num.sample}{num.sample is the number samples drawn to compute the sampling distibution.}
  \item{prop.sample}{prop.sample is the proportion of the original sample to include in the sampling distibution samples. Default is .9}
    \item{alpha}{The alpha value for confidence intervals. Default is .05.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{des}{Returns a data.frame containing the design matrix and additional columns for the predicted probabilities.}
  \item{boot.dist}{The full bootstrapped distribution for the probabilities.}
}
\references{
}
\author{
David Melamed
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("LF06travel")
require(Epi)
m1 <- clogistic(choice ~ train + bus + time  + invc, strata=id, data=LF06travel)
design <- data.frame(train=c(0,0,1),bus=c(0,1,0),time=200,invc=20)
design
margins.dat.clogit(m1,design)

ma1<-margins.dat.clogit(m1,design,run.boot="yes")
names(ma1)
ma1$des

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
