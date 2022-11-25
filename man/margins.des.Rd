\name{margins.des}
\alias{margins.des}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates a design matrix of idealized data for illustrating model predictions.
}
\description{
Create a data frame of idealized data for making model predictions/predicted margins.
Given a glm model object, a grid of indpendent variable values, and a list of any variables (factor variables in particular) to exclude from the design matrix, the function returns the design matrix as a data.frame object. All covariates are set to their means in the data used to estimate the model object. If there are factors in the model, they need to be excluded using the "excl" option. Note that, if the model is a multinomial regression model, as estimated via the nnet package, you need to also provide the data used in the nnet function that defined the model.
}
\usage{
margins.des(mod,ivs,excl="nonE",data)}
\arguments{
  \item{mod}{A model object. The model should be regression model for limited dependent variables, such as a logistic regression. Specifically, supported models include lm, glm, polr, multinom, vgam, zeroinf (pscl), hurdle (pscl) and zerotrunc (countreg). vgam is supported for partial proportional odds models, not models for count outcomes. zerotrunc is only supported with bootstrapped inference, and may take a while.
}
  \item{ivs}{This should be an 'expand.grid' statement including all desired variables and their corresponding levels in the design matrix.}
  \item{excl}{If you want to exclude covariates from the design matrix, you can list them here. This is designed to exclude factor variables from the design matrix, as they do not have means, but may be useful in other specialized cases. Default is "nonE," corresponding to excluding none of the variables.}
  \item{data}{If the model is a multinomial model, you also need to provide the data. This is because nnet objects do not include the relevant information for computing the means of covariates.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{design}{Returns a data.frame containing the design matrix for model predictions.}
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
data("Mize19AH")
m1 <- glm(alcB ~woman*parrole + age + race2 + race3 + race4 + income + ed1 + ed2 + ed3 + ed4,family="binomial",data=Mize19AH)
des1 <- margins.des(m1,expand.grid(parrole=1,woman=1))
des1
des2<-margins.des(m1,expand.grid(woman=c(0,1),parrole=c(0,1)))
des2
des3 <- margins.des(m1,expand.grid(age=seq(20,75,5),parrole=c(0,1)))
des3
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
