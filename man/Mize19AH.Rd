\name{Mize19AH}
\alias{Mize19AH}
\docType{data}
\title{
Add-Health Data analzed in Mize (2019)
}
\description{
Mize (2019) illustrates how to establish moderation in the context of regression models for limited dependent variables. He illustrates using AddHealth data and provides Stata code to replicate the results. Catregs functions can replicate these results in R.
}
\usage{data("Mize19AH")}
\format{
  A data frame with 4307 observations on the following 29 variables.
  \describe{
    \item{\code{AID}}{a numeric vector}
    \item{\code{race}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{educ}}{a numeric vector}
    \item{\code{degree}}{a numeric vector}
    \item{\code{college}}{a numeric vector}
    \item{\code{health}}{a numeric vector}
    \item{\code{role}}{a numeric vector}
    \item{\code{workrole}}{a numeric vector}
    \item{\code{parrole}}{a numeric vector}
    \item{\code{income}}{a numeric vector}
    \item{\code{wages}}{a numeric vector}
    \item{\code{logwages}}{a numeric vector}
    \item{\code{depB}}{a numeric vector}
    \item{\code{alcB}}{a numeric vector}
    \item{\code{woman}}{a numeric vector}
    \item{\code{edyrs}}{a numeric vector}
    \item{\code{whiteB}}{a numeric vector}
    \item{\code{X_est_prno}}{a numeric vector}
    \item{\code{X_est_prpar}}{a numeric vector}
    \item{\code{X_est_alcedmod}}{a numeric vector}
    \item{\code{X_est_alcmod}}{a numeric vector}
    \item{\code{race2}}{a numeric vector}
    \item{\code{race3}}{a numeric vector}
    \item{\code{race4}}{a numeric vector}
    \item{\code{ed1}}{a numeric vector}
    \item{\code{ed2}}{a numeric vector}
    \item{\code{ed3}}{a numeric vector}
    \item{\code{ed4}}{a numeric vector}
  }
}
\source{
Mize, Trenton D. 2019. "Best Practices for Estimating, Interpreting, and
Presenting Nonlinear Interaction Effects" Sociological Science 6: 81-117.
}
\examples{
data(Mize19AH)
head(Mize19AH)
}
\keyword{datasets}
