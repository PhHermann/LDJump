\name{getPhi}
\alias{getPhi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary Statistics to estimate recombination from PhiPack
%%  ~~function to do ... ~~
}
\description{
This functions calls the PhiPack software and extracts the four summary statistics MaxChi, NSS and the mean and the variance of Phi.
}
\usage{
getPhi(seqName, pathPhi, out, rm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{seqName}{
   A character string containing the full path and the name of the sequence file in \code{fasta} of \code{vcf} format. It is necessary to add the extension ("fileName.fa", "fileName.fasta", "fileName.vcf") in order to run \code{LDJump}. In case that \code{format} equals to \code{DNABin} the seqName equals to the name of the \code{DNABin}-object (without any extension).
}
  \item{pathPhi}{
  A character string containing the path to PhiPack. This path and the installation of \href{https://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz}{PhiPack} is necessary for the computation of the package.
  }
    \item{out}{
 an optional character string: by default an empty string "". Can be set to any user-defined string in order to rename all output files used within \code{LDJump} and \code{PhiPack}. This parameter enables to run \code{LDJump} from the same directory without creating interfering files in the working directory.
 }
 \item{rm}{
 an optional logical value: by default \code{TRUE} such that the internally produced fasta file as well as the output file are deleted shortly before finishing the function. This option is added in order to avoid deleting a file of interest when running the function \code{gethi} outside \code{LDJump}.
 }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector is returned containing the four summary statistics MaxChi, NSS and the mean and the variance of Phi.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bruen, T., Phillipe, H. and Bryant, D. 2006. A quick and robust statistical test to detect the presence of recombination. Genetics 172, 2665–2681.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LDJump}}, \code{\link{vcfR_to_fasta}}, \code{\link{get_smuce}}
}

\examples{
## The function is currently defined as
##getPhi(seqName = seqName, pathPhi = pathPhi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
