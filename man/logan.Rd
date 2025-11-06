\name{logan}
\alias{logan}
\docType{data}
\title{
Replication data for Logan's (1983) application of conditional logistic regression to mobility processes.
}
\description{
Replication data for Logan's (1983) application of conditional logistic regression to mobility processes.
}
\usage{data("logan")}
\format{
  A data frame with 4190 observations on the following 11 variables.
  \describe{
    \item{\code{occupation}}{respondent occupation with levels \code{farm} \code{operatives} \code{craftsmen} \code{sales} \code{professional}}
    \item{\code{focc}}{paternal occupation, i.e., father's occupation with levels \code{farm} \code{operatives} \code{craftsmen} \code{sales} \code{professional}}
    \item{\code{education}}{a numeric vector}
    \item{\code{race}}{a factor with levels \code{non-black} \code{black}}
    \item{\code{id}}{a numeric vector}
    \item{\code{tocc}}{a factor with levels \code{farm} \code{operatives} \code{craftsmen} \code{sales} \code{professional}}
    \item{\code{case}}{a numeric vector}
    \item{\code{craftsmen}}{a numeric vector}
    \item{\code{farm}}{a numeric vector}
    \item{\code{operatives}}{a numeric vector}
    \item{\code{professional}}{a numeric vector}
  }
}
\references{
Logan, John A. 1983. “A Multivariate Model for Mobility Tables.” American Journal of Sociology 89(2):324–349.
}
\examples{
data(logan)
head(logan)
}
\keyword{datasets}
