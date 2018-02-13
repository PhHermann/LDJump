\name{list.quantile.regs}
\alias{list.quantile.regs}
\docType{data}
\title{Quantile Regressions for Bias Correction
%%   ~~ data name/kind ... ~~
}
\description{
This data set contains a list of quantile regression models (\link[quantreg]{rq}) for bias correction. In total nine regression models are saved in a list form, where we recommend to use the 0.35 for the correction.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("list.quantile.regs")}
\format{
  A list object of length 9, containing quantiles in sequences between 0.1 and 0.5 with 0.05 distances.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(list.quantile.regs)
##### Do not run these examples                                           #####
##### In get_smuce the function is called as follows                      #####
##### pr1 = predict(mod.full,help); pr1[is.na(pr1)] = -1/gam;             #####
##### ind.q = which(seq(0.1, 0.5, by = 0.05) == quant)                    #####
##### pr.cor = predict(list.quantile.regs[[ind.q]], data.frame(x = pr1)   #####
}
\keyword{datasets}
\keyword{list}
\keyword{models, regression}
