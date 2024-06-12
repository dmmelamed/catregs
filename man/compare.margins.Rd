\name{compare.margins}
\alias{compare.margins}
\title{Compares two marginal effects (MEMs or AMEs). Estimate of uncertainty is from a simulated draw from a normal distribution.
}
\description{
Given two marginal effects (AMEs or MEMs), as estimated via the margins package or via first.diff.fitted, this function simulates draws from the distribution of MEs defined by the estimates and their standard error, and computes the overlap in the two distributions. The p-value refers to proportion of times the two draws overlapped.
}
\usage{
compare.margins(margins,margins.ses,seed=1234,rounded=3,nsim=10000)
}
\arguments{
  \item{margins}{The two marginal effects that you want to compare.}
  \item{margins.ses}{The standard errors for the marginal effects you want to compare.}
  \item{seed}{Random number seed so that results are reproducible.}
  \item{rounded}{The number of decimal places to round the output. The default is 3.}
  \item{nsim}{The number of simulated AMEs to draw from each distribution. The default is 10,000.}
}
\value{
  \item{differnce}{The observed difference in the two AMEs.}
  \item{p.value}{The p-value associated with the difference. This is the proportion of the simulated sample when the MEs overlapped.}}
\author{
David Melamed
}
\examples{
data("essUK")
m1 <- glm(safe ~ religious + minority*female + age,data=essUK,family="binomial")
des<-margins.des(m1,expand.grid(minority=c(0,1),female=c(0,1)))
des
ma1<-suppressWarnings(as.data.frame(marginaleffects::avg_slopes(m1,
variables="female",newdata=data.frame(minority=0,religious=3.6024,age=53.146))))
ma2<-as.data.frame(marginaleffects::avg_slopes(m1,variables="female",
newdata=data.frame(minority=1,religious=3.6024,age=53.146)))
cames <- rbind(ma2,ma1)
compare.margins(margins=cames$estimate,margins.ses=cames$std.error)
}
