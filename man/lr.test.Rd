\name{lr.test}
\alias{lr.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{LR test of nested models.
%%  ~~function to do ... ~~
}
\description{
Given two model objects estimated with a logLik() solution, the function computes an LR test. The order of models, in terms of full or reduced, does not matter.
}
\usage{
lr.test(full.model,reduced.model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full.model}{A model object with unconstrained parameters. The order does not matter however.}
  \item{reduced.model}{A model object with constrained parameters. The order does not matter however.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{LL Full}{The solution to the log-likelihood function of the full model.}
  \item{LL Reduced}{The solution to the log-likelihood function of the reduced model.}
  \item{G2/LR Statistic}{The test statistic. Computed as -2 x (ll{Full} - ll{Reduced}).}
  \item{DF}{The DF for the test. The number of parameters constrained to zero in the reduced model.}
  \item{p-value}{The p-value associated with the LR test. The null hypothesis is that the constrained parameters are jointly equal to zero.}
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
m0 <- glm(alcB ~woman + parrole + age + race2 + race3 + race4 + income + ed1 + ed2 + ed3 + ed4,family="binomial",data=Mize19AH)
lr.test(m1,m0)
lr.test(m0,m1) # Same Result
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
