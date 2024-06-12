\name{essUK}
\alias{essUK}
\docType{data}
\title{
A subset of data from the European Social Survey
}
\description{
These are data from respondents in the United Kingdom from the European Social Survey. They are used to illustrate regression models for limited dependent variables.
}
\usage{data("essUK")}
\format{
  A data frame with 2204 observations on the following 45 variables.
  \describe{
    \item{\code{country}}{a character vector}
    \item{\code{can.trust.people}}{a numeric vector}
    \item{\code{people.try.fair}}{a numeric vector}
    \item{\code{say.in.govt}}{a factor with levels \code{Not at all} \code{Very little} \code{Some} \code{A lot} \code{A great deal}}
    \item{\code{trust.legal.sys}}{a numeric vector}
    \item{\code{trust.police}}{a numeric vector}
    \item{\code{vote}}{a factor with levels \code{Yes} \code{No} \code{Not eligible to vote}}
    \item{\code{conservative}}{a numeric vector}
    \item{\code{life.satisfaction}}{a numeric vector}
    \item{\code{immigration.good.economy}}{a numeric vector}
    \item{\code{happy}}{a numeric vector}
    \item{\code{important.matters.people}}{a character vector}
    \item{\code{walk.alone.dark}}{a factor with levels \code{Very unsafe} \code{Unsafe} \code{Safe} \code{Very safe}}
    \item{\code{religious}}{a numeric vector}
    \item{\code{ethnic.minority}}{a character vector}
    \item{\code{num.children}}{a numeric vector}
    \item{\code{ideal.age.parent}}{a numeric vector}
    \item{\code{household.size}}{a numeric vector}
    \item{\code{gender}}{a character vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{marital}}{a character vector}
    \item{\code{education}}{a numeric vector}
    \item{\code{employment}}{a character vector}
    \item{\code{income.decile}}{a numeric vector}
    \item{\code{father.education}}{a character vector}
    \item{\code{father.education.num}}{a numeric vector}
    \item{\code{everyone.job.wanted}}{a numeric vector}
    \item{\code{income.fairness}}{a factor with levels \code{Low, extremely unfair} \code{Low, very unfair} \code{Low, somewhat unfair} \code{Low, slightly unfair} \code{Fair} \code{High, slightly unfair} \code{High, somewhat unfair} \code{High, very unfair} \code{High, extremely unfair}}
    \item{\code{under.over.paid}}{a factor with levels \code{Underpaid} \code{Right amount} \code{Overpaid}}
    \item{\code{income.fairness.num}}{a numeric vector}
    \item{\code{wealth.diff.fair}}{a factor with levels \code{Small, extremely unfair} \code{Small, very unfair} \code{Small, somewhat unfair} \code{Small, slightly unfair} \code{Fair} \code{Large, slightly unfair} \code{Large, somewhat unfair} \code{Large, very unfair} \code{Large, extremely unfair}}
    \item{\code{wealth.differences}}{a factor with levels \code{Too little} \code{Just right} \code{Too much}}
    \item{\code{gdp}}{a numeric vector}
    \item{\code{urban.population}}{a numeric vector}
    \item{\code{unemployment}}{a numeric vector}
    \item{\code{alcolhol.consumption}}{a numeric vector}
    \item{\code{suicide.rate}}{a numeric vector}
    \item{\code{safe}}{a numeric vector}
    \item{\code{minority}}{a numeric vector}
    \item{\code{female}}{a numeric vector}
    \item{\code{divorced}}{a numeric vector}
    \item{\code{married}}{a numeric vector}
    \item{\code{widow}}{a numeric vector}
    \item{\code{highinc}}{a numeric vector}
    \item{\code{age2}}{a numeric vector}
  }
}
\source{
European Social Survey.
}
\examples{
data(essUK)
head(essUK)}
\keyword{datasets}
