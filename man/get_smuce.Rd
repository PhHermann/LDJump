\name{get_smuce}
\alias{get_smuce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Segmentation Algorithm to Estimate Breakpoints in the Recombination Map
%%  ~~function to do ... ~~
}
\description{
First, the recombination rates per segment are computed based on the regression model (generalized additive models) as well as the bias correction.
    Consequently, we apply SMUCE (simultaneous multiscale change-point estimator) of Frick (2014) and Futschik et al. (2014) to estimate locations and breakpoints in the recombination map.
    Under a specific type-I error probability alpha the number of distinct segments with respect to the recombination rate is not overestimated.
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_smuce(help, segs, alpha, ll, quant = 0.35, rescale, constant, demography, regMod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{help}{
  a matrix containing a set of summary statistics is calculated in the function summary_statistics. These values are used in the regression model to calculate the (constant) recombination rates.
%%     ~~Describe \code{help} here~~
}
  \item{segs}{
  A (non-negative) integer which reflects the number of segments considered. It is calculated in the program based on the user-defined \code{segLength}.
%%     ~~Describe \code{segs} here~~
}
  \item{alpha}{
  A value from the interval (0,1) for the type-I error probability used in the segmentation algorithm. We recommend to use 0.05. We enabled to estimate the recombination map efficiently (without recalculating all summary statistics) under several type-I errors when \code{LDJump} is applied with a vector of type-I error probabilities.
%%     ~~Describe \code{alpha} here~~
}
  \item{ll}{
  A (non-negative) integer which reflects the total sequence length of the sequences under study.
%%     ~~Describe \code{ll} here~~
}
  \item{quant}{
  A value between 0.1 and 0.5 with 0.05 distances in between which reflects the quantile used in the quantile regression. We recommend to use the 0.35 quantile.
%%     ~~Describe \code{quant} here~~
}
  \item{rescale}{
  an optional logical value: if \code{TRUE}, it rescales the sequence length of the output of SMUCE to a range from 0 to 1.
%%     ~~Describe \code{rescale} here~~
}
  \item{constant}{
  an optional logical value: by default \code{FALSE} estimating variable recombination rates. If \code{TRUE}, the constant recombination rate for the full sequence is estimated.
  }
 \item{demography}{
  an optional character value: by default an empty string ("") indicates that the recombination rate estimation is estimated under neutrality. If \code{"b"} the regression model estimated based on samples from populations under a bottleneck is used. If \code{"g"} the regression model estimated based on samples from populations under population growth is used. If \code{"d"}, the regression model estimated based on samples from populations under demography (combination of samples of under growth and bottleneck) is used.
}
 \item{regMod}{
 an optional character string: for the default empty string "" *LDJump* uses an existing regression model (constant population size or simple demography example, depending on \code{demography}). In oder to use the regression model estimated by user input demography, then this variable should equal to the name of the regression object. Please see the examples for more details.
 }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{seq.full.cor}{The final estimate of the recombination map. Depiction with plot-function of stepR package. }
  \item{pr.full.cor}{A vector of (constant) estimates of the recombination rate per segment. }

%% ...
}
\references{
Frick, K., Munk, A., and Sieling, H. (2014). Multiscale change-point inference.
Journal of the Royal Statistical Society: Series B, 76(3), 495–580.

Futschik, A., Hotz, T., Munk, A., and Sieling, H. (2014). Multiscale DNA
partitioning: Statistical evidence for segments. Bioinformatics, 30(16), 2255–2262.

Hermann, P., Heissl, A., Tiemann-Boege, I., and Futschik, A. (2019), LDJump:
Estimating Variable Recombination Rates from Population Genetic Data.
Mol Ecol Resour. \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12994?af=R}{doi:10.1111/1755-0998.12994}.

}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LDJump}}, \code{\link{vcfR_to_fasta}}, \code{\link{getPhi}}, \code{\link{summary_statistics}}, \code{\link[stepR]{stepFit}}, \code{\link[quantreg]{rq}}, \code{\link[mgcv]{gam}}
}

\examples{
##### Do not run these examples                           #####
##### In LDJump.R the function is called as follows       #####
##### get_smuce(help, segs, alpha,ll,list.quantile.regs)  #####

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{htest}
