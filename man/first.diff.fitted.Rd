\name{first.diff.fitted}
\alias{first.diff.fitted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the first difference in fitted values. Inference in supported via the delta method or bootstrapping.
}
\description{
first.diff.fitted computes the first difference between fitted values from a regression model.
}
\usage{
first.diff.fitted(mod,design.matrix,compare,alpha=.05,rounded=3,bootstrap="no",num.sample=1000,prop.sample=.9,data)
}
\arguments{
  \item{mod}{A model object. The model should be regression model for limited dependent variables, such as a logistic regression.}
  \item{design.matrix}{Design matrix of values for the independent variables in the regression model.}
  \item{compare}{A pair of rows in the design matrix to use for computing the fitted values. The first difference between the fitted values is then computed. For example, compare=c(4,2) means to compute the difference in the fitted values between predictions for row 4 of the design matrix and row 2 of the design matrix. Only two rows may be compared at a time.}
  \item{alpha}{The alpha value for confidence intervals. Default is .05.}
  \item{rounded}{The number of decimal places to round the output. The default is 3.}
  \item{bootstrap}{By default, inference is based on the Delta Method, as implemented in the marginaleffects package. Alternatively, inference can be based upon a bootstrapped sampling distirbution. To do so, change this to "yes"}
    \item{num.sample}{num.sample is the number samples drawn to compute the sampling distibution when using bootstrapping. Default is 1,000}
  \item{prop.sample}{prop.sample is the proportion of the original sample to include in the sampling distibution samples when using bootstrapping. Default is .9}
  \item{data}{For nonparametric inference, provide the data used in the original model statement.}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{out}{If using parametric inference (delta method): output is a dataframe including the difference in fitted values ("est"), the standard error ("std.error"), the lower limit ("ll"), and upper limit ("ul") of the confidence interval. Of course, ll and ul are based on the alpha level.
  If using nonparametric inference (bootstrapping): output is a list of objects. obs.diff is the observed difference in the response or fitted values. boot.dist is the sorted bootstrapped distribution of differences in the samples. mean.boot.dist is the average of the differences in the responses or fitted values. sd.boot.dist is the standard deviation of the sampling distribution. ci.95 is the Lower and Upper limits of the confidence interval; despite it's name, the confidence interval is based upon the alpha level. model.class is just the class of the model that was used to generate the fitted values.}
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
data("Mize19AH")
m1 <- glm(alcB ~woman*parrole + age + race2 + race3 + race4 + income + ed1 + ed2 + ed3 + ed4,family="binomial",data=Mize19AH)
des2<-margins.des(m1,expand.grid(woman=c(0,1),parrole=c(0,1)))
des2
first.diff.fitted(m1,des2,compare=c(4,2)) # Pr(Drink | Mothers) - Pr(Drink | Childless Women)
first.diff.fitted(m1,des2,compare=c(3,1)) # Pr(Drink | Fathers) - Pr(Drink | Childless Men)
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
