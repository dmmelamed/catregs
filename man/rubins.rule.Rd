\name{rubins.rule}
\alias{rubins.rule}
\title{Aggregate Standard Errors using Rubin's Rule.
}
\description{
The function takes a vector of standard error estimates and it pools them using Rubin's rule.
}
\usage{
rubins.rule(std.errors)
}
\arguments{
  \item{std.errors}{A vector of standard errors to be aggregated using Rubin's rule.}
}
\value{
  \item{r.r.std.error}{The aggregated standard error.}
}
\references{Rubin, Donald B. 2004. Multiple Imputation for Nonresponse in Surveys. Vol. 81. John Wiley & Sons.

}
\author{
David Melamed
}
