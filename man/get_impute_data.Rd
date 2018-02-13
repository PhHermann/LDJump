\name{get_impute_data}
\alias{get_impute_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data management to obtain data used for imputation
%%  ~~function to do ... ~~
}
\description{
This function obtains the neighbouring values for every segment without SNPs. The function either merges data of two or four neighbouring values. It is only used in the \code{impute} function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_impute_data(index, data,two=T,segs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{index}{
  this is a vector containing the integer number of the segments without SNPs.
%%     ~~Describe \code{x} here~~
}
  \item{data}{
  A data vector containing the estimated recombination rates per segment.
  }
  \item{two}{
  A logical parameter indicating whether two neighbouring values (if \code{TRUE}) or four neighbouring values (if \code{FALSE}) should be provided.
  }
    \item{segs}{
  A (non-negative) integer which reflects the number of segments considered. It is calculated in the program based on the user-defined \code{segLength}.
%%     ~~Describe \code{segs} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 The function returns a vector containing the corresponding data which can be used for imputation. Note that the ordering is made that values of distance 1 are listed first in case that four neighbouring values are returned.
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
\code{\link{impute}}
}
\examples{
##### Do not run these examples                                  #####
##### This command shows how it is used in the impute function   #####
##### sapply(index, get.impute.data, data, two = two)            #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{datasets}
