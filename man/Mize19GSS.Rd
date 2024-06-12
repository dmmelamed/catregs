\name{Mize19GSS}
\alias{Mize19GSS}
\docType{data}
\title{
General Social Survey Data analzed in Mize (2019)}
\description{
Mize (2019) illustrates how to establish nonlinear moderation in the context of regression models. He illustrates using General Social Survey (GSS) data and provides Stata code to replicate the results. Catregs functions can replicate these results in R.}
\usage{data("Mize19GSS")}
\format{
  A data frame with 19337 observations on the following 42 variables.
  \describe{
    \item{\code{nosameB}}{a numeric vector}
    \item{\code{sameokB}}{a numeric vector}
    \item{\code{polviews}}{a character vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{age10}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{id}}{a numeric vector}
    \item{\code{degree}}{a numeric vector}
    \item{\code{race}}{a numeric vector}
    \item{\code{partyid}}{a character vector}
    \item{\code{natspac}}{a character vector}
    \item{\code{natenvir}}{a character vector}
    \item{\code{natheal}}{a character vector}
    \item{\code{natcity}}{a character vector}
    \item{\code{natcrime}}{a character vector}
    \item{\code{natdrug}}{a character vector}
    \item{\code{nateduc}}{a character vector}
    \item{\code{natrace}}{a character vector}
    \item{\code{natarms}}{a character vector}
    \item{\code{nataid}}{a character vector}
    \item{\code{natfare}}{a character vector}
    \item{\code{health}}{a character vector}
    \item{\code{helpnotB}}{a character vector}
    \item{\code{conserv}}{a character vector}
    \item{\code{polviews3}}{a character vector}
    \item{\code{employed}}{a numeric vector}
    \item{\code{male}}{a numeric vector}
    \item{\code{woman}}{a numeric vector}
    \item{\code{white}}{a numeric vector}
    \item{\code{college}}{a numeric vector}
    \item{\code{married}}{a numeric vector}
    \item{\code{parent}}{a character vector}
    \item{\code{edyrs}}{a numeric vector}
    \item{\code{income}}{a numeric vector}
    \item{\code{hrswork}}{a character vector}
    \item{\code{parttime}}{a character vector}
    \item{\code{wages}}{a numeric vector}
    \item{\code{conviewSS}}{a numeric vector}
    \item{\code{year2}}{a numeric vector}
    \item{\code{yearcat}}{a numeric vector}
    \item{\code{year1976}}{a numeric vector}
    \item{\code{year1976.2}}{a numeric vector}
  }
}
\source{
Mize, Trenton D. 2019. "Best Practices for Estimating, Interpreting, and
Presenting Nonlinear Interaction Effects" Sociological Science 6: 81-117.
}
\examples{
data(Mize19GSS)
head(Mize19GSS)
}
\keyword{datasets}
