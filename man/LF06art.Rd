\name{LF06art}
\alias{LF06art}
\docType{data}
\title{
Data to replicate Long and Freese's (2006) count models (pp354-414)
}
\description{
For replication purposes between Stata and R. Long and Freese (2006) analyze these data to illustrate regression models for count dependent variables. }
\usage{data("LF06art")}
\format{
  A data frame with 915 observations on the following 6 variables.
  \describe{
    \item{\code{art}}{count response}
    \item{\code{fem}}{dummy for sex}
    \item{\code{mar}}{dummy for married}
    \item{\code{kid5}}{number of children under five}
    \item{\code{phd}}{a numeric vector}
    \item{\code{ment}}{a numeric vector}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Long, Scott J. and Jeremy Freese. 2006. "Regression Models for Categorical Dependent Variables Using Stata." Austin, TX: Stata Press
}

\examples{
data(LF06art)
head(LF06art)
}
\keyword{datasets}
