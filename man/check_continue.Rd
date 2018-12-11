\name{check_continue}
\alias{check_continue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checks whether there are SNPs in each segment
%%  ~~function to do ... ~~
}
\description{This function calculates the number of SNPs per segment. In case that there exist segments with less than 2 SNPs the user is asked for input to continue ("y") or not ("n"). In case that the user wants to continue, the recombination rates for segments without SNPs are estimated via imputation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
check_continue(seqName, segs, accept, format)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqName}{
%%     ~~Describe \code{seqName} here~~
     A character string containing the full path and the name of the sequence file in \code{fasta} of \code{vcf} format. It is necessary to add the extension ("fileName.fa", "fileName.fasta", "fileName.vcf") in order to run \code{LDJump}. In case that \code{format} equals to \code{DNABin} the seqName equals to the name of the \code{DNABin}-object (without any extension).
}
  \item{segs}{
  A (non-negative) integer which reflects the number of segments considered. It is calculated in the program based on the user-defined segment length.
  %%     ~~Describe \code{segs} here~~
}
  \item{accept}{
  an optional logical value: by default \code{FALSE} and \code{LDJump} checks for segments with less than 2 SNPs and requires user input to proceed. If set to \code{TRUE}, the user accepts that the rates for these segments (<= 1 SNP) are estimated via imputation.
  }
  \item{format}{
  A character string which can be \code{fasta}, \code{vcf}, or \code{DNAStringSet}. If \code{fasta} is used, the package will proceed with the computation of the recombination map. If \code{vcf}, the package will convert the data in \code{vcf} format to \code{fasta} format with the function \code{vcfR_to_fasta} and then proceed as in case \code{fasta}. For the last format the \code{seqName} must equal to the \code{DNABin}-object which contains the sequences.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 This function returns \code{TRUE} in case that all segments contain SNPs. It will also return \code{TRUE} if the user agrees to continue although there exist segments without SNPs. It returns \code{FALSE} if the user denies to continue due to segments without SNPs.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Paradis E., Claude J. & Strimmer K. 2004. APE:  analyses of phylogenetics and evolution
in R language. Bioinformatics 20: 289-290.
}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[ape]{read.FASTA}}, \code{\link[ape]{seg.sites}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##### Do not run these examples                                  #####
##### check_continue(seqName, segs = segs)                       #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
