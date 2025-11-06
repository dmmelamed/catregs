\name{LF06travel}
\alias{LF06travel}
\docType{data}
\title{
Travel time example data for alternative-specific outcomes.
}
\description{
Example data, also used in Long and Freese (2006), to illustrate conditional or fixed effects logistic regression. Also refered to as alternative-specific outcomes.
}
\usage{data("LF06travel")}
\format{
  A data frame with 456 observations on the following 13 variables.
  \describe{
    \item{\code{id}}{a numeric vector denoting nested units (individuals) or strata}
    \item{\code{mode}}{a numeric vector denoting mode of transit}
    \item{\code{train}}{a dummy variable for selecting the train}
    \item{\code{bus}}{a dummy variable for selecting the bus}
    \item{\code{car}}{a dummy variable for selecting a car}
    \item{\code{time}}{a numeric vector denoting transit time}
    \item{\code{invc}}{a numeric vector denoting invertng cost}
    \item{\code{choice}}{a numeric vector denoting the choice of travel, i.e. the dependent variable}
    \item{\code{ttme}}{a numeric vector}
    \item{\code{invt}}{a numeric vector}
    \item{\code{gc}}{a numeric vector}
    \item{\code{hinc}}{a numeric vector}
    \item{\code{psize}}{a numeric vector}
  }
}
\source{
Long, Scott J. and Jeremy Freese. 2006. "Regression Models for Categorical Dependent Variables Using Stata." Austin, TX: Stata Press
}
\examples{
data(LF06travel)
head(LF06travel)
}
\keyword{datasets}
