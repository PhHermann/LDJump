\name{mod.full}
\alias{mod.full}
\docType{data}
\title{Regression model to Estimate (constant) Recombination Rates from Population Genetic Summary Statistics
%%   ~~ data name/kind ... ~~
}
\description{
This data set contains a generalized additive regression model (\link[mgcv]{gam}) which estimates the constant recombination rate for a set of segments based on summary statistics.
%%  ~~ A concise  (1-5 lines) description of the dataset. ~~
}
\usage{data("mod.full")}
\format{
  A list containing all information on the regression model such as the coefficients, residuals, among others as usually for generalized additive models (\link[mgcv]{gam}) saved as an \code{R object}.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
  Wood, S.N. (2011) Fast stable restricted maximum   likelihood and marginal
  likelihood estimation of   semiparametric generalized linear models.
  Journal of the Royal Statistical Society (B) 73(1):3-36
}
\examples{
data(mod.full)
##### Do not run these examples                         #####
##### In get_smuce the function is called as follows    #####
##### pr1 = predict(mod.full,help)                      #####
}
\keyword{datasets}
\keyword{models, regression}
