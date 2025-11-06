\name{gss2016}
\alias{gss2016}
\docType{data}
\title{
Data from the 2016 General Social Survey.
}
\description{
Limited date from the 2016 General Social Survey on respondent and paternal class and occupational classifications.
}
\usage{data("gss2016")}
\format{
  A data frame with 12498 observations on the following 13 variables.
  \describe{
    \item{\code{pclass}}{a factor with levels \code{Unskilled Manual} \code{Skilled Manual} \code{Self-Employed} \code{Non-Manual/Service} \code{Professional, Lower} \code{Professional, Higher}}
    \item{\code{sclass}}{a factor with levels \code{Unskilled Manual} \code{Skilled Manual} \code{Self-Employed} \code{Non-Manual/Service} \code{Professional, Lower} \code{Professional, Higher}}
    \item{\code{educ}}{a numeric vector}
    \item{\code{race}}{a character vector}
    \item{\code{id}}{a numeric vector}
    \item{\code{occ2}}{a character vector}
    \item{\code{occ}}{a numeric vector}
    \item{\code{unskmanual}}{a numeric vector}
    \item{\code{skmanual}}{a numeric vector}
    \item{\code{selfemp}}{a numeric vector}
    \item{\code{service}}{a numeric vector}
    \item{\code{proflow}}{a numeric vector}
    \item{\code{profhigh}}{a numeric vector}
  }
}
\source{
The General Social Survey.
}
\examples{
data(gss2016)
head(gss2016)
}
\keyword{datasets}
