\name{diagn}
\alias{diagn}
\title{Computes diagnostics for generalized linear models.
}
\description{
Given a glm object, diagn returns case-level diagnostics.
For logistic, probit, Poisson, and negative binomial models, it returns Pearson residuals, standardized Pearson residuals, the diagonal of the hat matrix, delta-beta (Cook's D), and deviance residuals. For zero-inflated and hurdle models, it returns the Pearson residual and the observation number.
}
\usage{
diagn(model)
}
\arguments{
  \item{model}{A model object. The model should be regression model for limited dependent variables, such as a logistic regression.}
}
\value{
  \item{out}{The output is a dataframe of diagnostic statistics.
  For logit, probit, Poisson, and negative binomial models, the output includes the Pearson residual (pearsonres), the diagonal of the Hat matrix (h), the standardized Pearson residual (stdpres), the delta-beta statistic (deltabeta), the observation number (obs), and the deviance residual (devres).
  For zero-inflated and hurdle models, the output includes the Pearson residual (pearsonres), and the observation number (obs).}
}
\author{
David Melamed
}
\examples{
data("Mize19AH")
m1 <- glm(alcB ~woman*parrole + age + race2 + race3 +
race4 + income + ed1 + ed2 + ed3 + ed4,
family="binomial",data=Mize19AH)
head(diagn(m1))
}
