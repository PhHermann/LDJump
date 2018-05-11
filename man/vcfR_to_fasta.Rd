\name{vcfR_to_fasta}
\alias{vcfR_to_fasta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conversion of vcf to fasta Format
%%  ~~function to do ... ~~
}
\description{
This function enables to read vcfR files and convert them to necessary fasta files.
    Therefore, we recommend to provide a reference sequence from e.g. genome browser and the starting position.
    The default parameters are those of the vcfR package.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vcfR_to_fasta(seqName, refName = NULL, ext.ind = T, cons = F,
              ext.haps = T, start = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{seqName}{
%%     ~~Describe \code{seqName} here~~
   A character string containing the full path and the name of the sequence file in vcf format.
%%     ~~Describe \code{seqName} here~~
}
  \item{refName}{
  An (optional) path to the reference sequence for the region of interest downloaded from e.g. http://phase3browser.1000genomes.org/index.html. Only to be used in case that \code{format == "vcf"}.
  }
%%     ~~Describe \code{refName} here~~
  \item{ext.ind}{
  See package \pkg{vcfR} for details (\code{\link[vcfR]{vcfR2DNAbin}}, \code{extract.indels})
%%     ~~Describe \code{ext.ind} here~~
}
  \item{cons}{
  See package \pkg{vcfR} for details (\code{\link[vcfR]{vcfR2DNAbin}}, \code{consensus})
%%     ~~Describe \code{cons} here~~
}
  \item{ext.haps}{
  See package \pkg{vcfR} for details (\code{\link[vcfR]{vcfR2DNAbin}}, \code{extract.haps})
%%     ~~Describe \code{ext.haps} here~~
}
  \item{start}{
  An (optional) integer value which reflects the starting position of the sequences in bp. Only to be used in case that \code{format == "vcf"}.
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
 A print command provides information that the file is converted.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Knaus BJ and Grünwald NJ (2017). VCFR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources, 17(1), pp. 44-53. ISSN 757, <URL:
http://dx.doi.org/10.1111/1755-0998.12549>.
}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LDJump}}, \code{\link{summary_statistics}}, \code{\link{fgt_rrate_dpr}}, \code{\link{get_smuce}}, \code{\link[vcfR]{vcfR2DNAbin}}
}

\examples{
##### Do not run these examples                                #####
##### vcfR_to_fasta (seqName, refName, ext.ind = T, cons = F,  #####
#####                ext.haps = T, start = 1)                  #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
