\name{fgt_rrate_dpr}
\alias{fgt_rrate_dpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Four Gametes Test, R^2 and LD' Calculation
%%  ~~function to do ... ~~
}
\description{
This helper function calculates three summary statistics of the regression model. Here, the four gametes test, R^2, and LD' are calculated for each pair of sites and returned to its calling function (summary_statistics)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fgt_rrate_dpr(x, y, data1, data2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Site 1
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  Site 2
%%     ~~Describe \code{y} here~~
}
  \item{data1}{
  Data set to calculate the first two summary statistics
%%     ~~Describe \code{data1} here~~
}
  \item{data2}{
  Data set to calculate the four gametes test
%%     ~~Describe \code{data2} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  A vector is returned containing three values for:
  \item{fgt}{An indicator value whether the four gametes test indicates a recombination event}
  \item{R^2}{for the pair of sites x and y using \code{\link[genetics]{diseq}}, \code{\link[genetics]{genotype}}}
  \item{LD'}{for the pair of sites x and y using \code{\link[genetics]{diseq}}, \code{\link[genetics]{genotype}}}
%% ...
}

\references{
Gregory Warnes, with contributions from Gregor Gorjanc, Friedrich Leisch and Michael Man. (2013).
genetics: Population Genetics. R package version 1.3.8.1.
}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LDJump}}, \code{\link{vcfR_to_fasta}}, \code{\link{summary_statistics}}, \code{\link{get_smuce}}, \code{\link[genetics]{diseq}}, \code{\link[genetics]{genotype}}
}

\examples{
##### Do not run these examples                                          #####
##### In summary_statistics.R the function is called as follows          #####
##### helper = mapply(fgt_rrate_dpr, x = indices[,1], y = indices[,2],   #####
#####                 MoreArgs = list(data1 = g2dftemp, data2 = samp))   #####

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

