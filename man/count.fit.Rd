\name{count.fit}
\alias{count.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits four different count models and compares them.
}
\description{
Given a Poisson model object, count.fit fits Poisson, negative binomial, zero-inflated Poisson, and zero-inflated negative binomial models to the data. It reports results of Vuong tests between the zero-inflated and non-zer-inflated models, summarizes the information criteria of the four models, summarizes the model output of the four models, creates a ggplot object of coefficient plots for each model, and creates a ggplot object of model residuals.
}
\usage{
count.fit(m1,y.range,rounded=3,use.color="yes")
}
\arguments{
  \item{m1}{A Poisson regression model, as estimated via the glm function.}
  \item{y.range}{The observed response range for the count outcome. For example, if the observed range is 0 to 18, this would be 0:18}
    \item{rounded}{The number of decimal places to round the output. The default is 3.}
    \item{use.color}{Whether to use color in the ggplot objects. Default is "yes"}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ic}{A data.frame summarizing the information criteria for the four models. Bayesian and Akaike's informaiton criteria are included.}
  \item{models}{A summary of the model estimates, including coefficients and standard errors.}
  \item{pic}{A ggplot object illustrating model residuals for each type of model.}
  \item{models.pic}{A ggplot object of coefficient plots from each type of model.}
}
\references{
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
data("LF06art")
p1 <- glm(art ~ fem + mar + kid5 + phd + ment , family = "poisson", data = LF06art)
table(LF06art$art)
fit<-count.fit(p1,0:19)
names(fit)
fit$ic
fit$pic
fit$models
fit$models.pic
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
