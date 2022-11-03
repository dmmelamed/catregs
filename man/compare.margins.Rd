\name{compare.margins}
\alias{compare.margins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compares two marginal effects (MEMs or AMEs). Estimate of uncertainty is from a simulated draw from a normal distribution.
%%  ~~function to do ... ~~
}
\description{
Given two marginal effects (AMEs or MEMs), as estimated via the margins package or via first.diff.fitted, this function simulates draws from the distribution of MEs defined by the estimates and their standard error, and computes the overlap in the two distributions. The p-value refers to proportion of times the two draws overlapped.
}
\usage{
compare.margins(margins,margins.ses,seed=1234,rounded=3,nsim=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{margins}{The two marginal effects that you want to compare.}
  \item{margins.ses}{The standard errors for the marginal effects you want to compare.}
  \item{seed}{Random number seed so that results are reproducible.}
  \item{rounded}{The number of decimal places to round the output. The default is 3.}
  \item{nsim}{The number of simulated AMEs to draw from each distribution. The default is 10,000.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{differnce}{The observed difference in the two AMEs.}
  \item{p.value}{The p-value associated with the difference. This is the proportion of the simulated sample when the MEs overlapped.}}
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
data("essUK")
m1 <- glm(safe ~ religious + minority*female + age,data=essUK,family="binomial")
require(margins)
ma1 <- summary(margins(m1,variables="female",at=list(minority=0)))
ma2 <- summary(margins(m1,variables="female",at=list(minority=1)))
cames <- rbind(ma2,ma1)
compare.margins(margins=cames$AME,margins.ses=cames$SE)
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
