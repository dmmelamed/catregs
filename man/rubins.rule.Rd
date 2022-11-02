\name{rubins.rule}
\alias{rubins.rule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Aggregate Standard Errors using Rubin's Rule.
%%  ~~function to do ... ~~
}
\description{
The function takes a vector of standard error estimates and it pools them using Rubin's rule.
}
\usage{
rubins.rule(std.errors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{std.errors}{A vector of standard errors to be aggregated using Rubin's rule.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{r.r.std.error}{The aggregated standard error.}
}
\references{Rubin, Donald B. 2004. Multiple Imputation for Nonresponse in Surveys. Vol. 81. John Wiley & Sons.

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
See Chapter 11 of the book.
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
