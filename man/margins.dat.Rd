\name{margins.dat}
\alias{margins.dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add model predictions, standard errors and confidence intervals to a design matrix for a model object.
}
\description{
Given a model object and a design matrix, this creates a data.frame of the design matrix, with model predictions, standard errors and lower/upper limits of confidence intervals around the predictions. This is a wrap around function for calls to emmeans; it adjusts the emmeans equation to return fitted values on the response scale.

Supported models include OLS regression via lm, logistic regression via glm, Poisson regression via glm, negative binomial regression via MASS:glm.nb, ordinal logistic regression via MASS::polr, partial proportional odds models via vgam::vglm, multinomial logistic regression via nnet::multinom, zero-inflated Poisson or negative binomial regression via pscl::zeroinfl, zero-truncated Poisson or negative binomial regression via countreg::zerotrunc, hurdle Poisson or negative binomial regression via pscl::hurdle, linear mixed effects models with or without serial correlation via nlme::lme, linear mixed effects models via lme4/lmerTest::lmer, mixed effects logistic regression via lme4/lmerTest::glmer, mixed effects Poisson regression via lme4/lmerTest::glmer, mixed effects negative binomial regression via lme4/lmerTest::glmer.nb, and mixed effects ordinal logistic regression via ordinal::clmm. zerotrunc is only supported with bootstrapped inference, and may take a while as a result.

Given one of these model objects and an appropiate design matrix, the function detects the model response type and generates fitted values on the response scale. For example, a logistic regression model returns predicted probabilities, and a Poisson model returns the fitted counts. In addition to the fitted values, the function returns the delta method standard error for the fitted value and a confidence interval. The confidence interval is 95 percent by default, but that may be changed by the user.
}
\usage{
margins.dat(mod,des,alpha=.05,rounded=3,cumulate="no",pscl.data=data,num.sample=1000,prop.sample=.9,seed=1234)}
\arguments{
  \item{mod}{A regression model object.}
  \item{des}{Design matrix of values for the independent variables in the regression model.}
  \item{alpha}{The alpha value for confidence intervals. Default is .05.}
  \item{rounded}{The number of decimal places to round the output. The default is 3.}
  \item{cumulate}{Whether the fitted values should reflect cumulative probabilities. Default is "no." Intended for predicted probabilities drawn from ordinal logistic regression models (polr), ordinal logistic regression models assuming partial proportional odds (vgam), or mixed effects ordinal logistic regression (clmm).}
  \item{pscl.data}{If generating predicted counts from Zero-Inflated models (either Poisson or negative binomial), you need to include the data that was specified in the model statement, i.e., the data in the "mod" object.}
  \item{num.sample}{Inference for zero truncated models, as estimated by the countreg package, is only supported via non-parametric/bootstrap inference. num.sample is the number samples drawn to compute the sampling distibution.}
  \item{prop.sample}{Inference for zero truncated models, as estimated by the countreg package, is only supported via non-parametric/bootstrap inference. prop.sample is the proportion of the original sample to include in the sampling distibution samples. Default is .9}
  \item{seed}{For models using bootstrapped inference. The seed ensures reproducible results across runs. Default is 1234, but may be changed.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{marginsdat}{Returns a data.frame containing the design matrix and additional columns for the fitted value on the response scale, the delta method standard error (except zero-truncated models, which are bootstrapped), and the lower/upper limits on confidence intervals around the fitted value.}
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
margins.dat(m1,des2,rounded=5)
des1 <- margins.des(m1,expand.grid(parrole=1,woman=1))
margins.dat(m1,des1,rounded=5)
des3 <- margins.des(m1,expand.grid(age=seq(20,75,5),parrole=c(0,1)))
a<- margins.dat(m1,des3,rounded=5)
a
ggplot(a,aes(age,y=fitted,ymin=ll,ymax=ul,group=as.factor(parrole),color=as.factor(parrole),linetype=as.factor(parrole))) +
  theme_bw() + geom_line() + geom_ribbon(alpha=.2) +
  theme(legend.position="bottom") +labs(x="Age",y="Pr(Alcohol Use)",color="Parent",linetype="Parent")
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
