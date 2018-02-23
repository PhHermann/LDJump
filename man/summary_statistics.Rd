\name{summary_statistics}
\alias{summary_statistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary Statistics per Segment
%%  ~~function to do ... ~~
}
\description{
This function computes summary statistics for every segment of the sequence.
    Sequence files are generated within this function which are then used by \href{https://github.com/auton1/LDhat}{LDhat} and other packages to estimate all necessary parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary_statistics(x, s, segLength, segs, seqName, nn, thth,
                   cor = 1, pathLDhat, status, polyThres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
An integer control variable for the considered segment of the DNA sequence.
}
  \item{s}{
%%     ~~Describe \code{s} here~~
An \code{\link[Biostrings]{XStringSet}} object which is read by \code{\link[Biostrings]{readDNAStringSet}}
}
  \item{segLength}{
%%     ~~Describe \code{segLength} here~~
   An integer value for the length of the segments, provided by the user. The default value of 1000 is our recommended value (1kb). The number of resulting segments, based on the sequence length is calculated within the funtion.
}
  \item{segs}{
  A (non-negative) integer which reflects the number of segments considered. It is calculated in the program based on the user-defined \code{segLength}.
%%     ~~Describe \code{segs} here~~
}
\item{seqName}{
%%     ~~Describe \code{seqName} here~~
   A character string containing the path and the name of the sequence file in fasta of vcf format.
%%     ~~Describe \code{seqName} here~~
}
  \item{nn}{
An integer which reflects the number of individuals (more precisely sequences) of the population to be analyzed. In case of diploid samples this is twice the number of individuals.
%%     ~~Describe \code{nn} here~~
}
  \item{thth}{
  A numeric value for \eqn{\theta} used in the lookup tables of \href{https://github.com/auton1/LDhat}{LDhat}.
%%     ~~Describe \code{thth} here~~
}
  \item{cor}{
  An integer value which reflects the number of cores on which \href{https://github.com/auton1/LDhat}{LDhat} should be run. We recommend to keep here 1 core.
  }
  \item{pathLDhat}{
  A character string containing the path to LDhat. This path and the installation of \href{https://github.com/auton1/LDhat}{LDhat} is necessary for the computation of the package.
  }
  \item{status}{
  an optional logical value: by default \code{TRUE} such that the current processing status of the segments is printed.
  }
  \item{polyThres}{
  a numeric value between 0 and 1. Used in data manipulation function \code{DNAbin2genind}: the minimum frequency of a minor allele for a locus to be considered as polymorphic (default to 0).
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
 This function returns a concatenated vector of two separately used vectors (scalars) of summary statistis as:
  \item{stats}{A vector of summary statistics. Returned with the value of hahe.}
  \item{hahe}{The haplotype heterozygosity of the considered segment. Returned with stats.}
  \item{tajd}{Tajima's D. Only used in the regression model for demography.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Auton, A. and McVean, G. (2007). Recombination rate estimation in the presence
of hotspots. Genome Research, 17(8), 1219–1227.

Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new   tools for the analysis of
genome-wide SNP data. Bioinformatics. doi: 10.1093/bioinformatics/btr521

McVean, G. A. T., Myers, S. R., Hunt, S., Deloukas, P., Bentley, D. R., and Donnelly,
P. (2004). The fine-scale structure of recombination rate variation in the human
genome. Science, 304(5670), 581–584.

Paradis E., Claude J. & Strimmer K. 2004. APE:  analyses of phylogenetics and evolution
in R language. Bioinformatics 20: 289-290.


}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LDJump}}, \code{\link{vcfR_to_fasta}}, \code{\link{fgt_rrate_dpr}}, \code{\link{get_smuce}}, \code{\link[Biostrings]{readDNAStringSet}}, \code{\link[adegenet]{DNAbin2genind}}
}

\examples{
##### Do not run these examples                                                      #####
##### In LDJump.R the function is called as follows                                  #####
##### sapply(1:segs,summary_statistics,s=s,segs=segs,seqName=seqName,nn=nn,ll = ll)  #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")

