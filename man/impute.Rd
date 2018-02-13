\name{impute}
\alias{impute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Imputation of estimated recombination rates for segments without SNPs
%%  ~~function to do ... ~~
}
\description{
This function recursively imputes the recombination rate for missing segments. First it imputes the mean of the two neighbouring segments. In case that one of these segments is also missing, it then imputes the weighted mean of the four neighbouring segments putting higher weights to the closer segments. Exceptions were made for e.g. the first and last segment in the sequence. It also imputes those positions first, where more information is already available.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
impute(data, index, two, segs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
  A data vector containing the estimated recombination rates per segment.
  }
  \item{index}{
  this is a vector containing the integer number of the segments without SNPs.
%%     ~~Describe \code{x} here~~
  }
  \item{two}{
  A logical parameter indicating whether two neighbouring values (if \code{TRUE}) or four neighbouring values (if \code{FALSE}) should be provided.
  }
  \item{segs}{
  A (non-negative) integer which reflects the number of segments considered. It is calculated in the program based on the user-defined \code{segLength}.
%%     ~~Describe \code{segs} here~~
}
}
\details{
The function calls itself after every imputation step trying to impute based on two neighbouring segments.
}
\value{
\item{data}{This vector contains the estimated recombination rates including the imputed values for the segments without SNPs.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{get_impute_data}}
}
\examples{
##### Do not run these examples                                   #####
##### This command shows how it is used in the get_smuce function #####
##### pr.cor.nat = impute(pr.cor.nat, ind, two = T)               #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
\keyword{datasets}% __ONLY ONE__ keyword per line
