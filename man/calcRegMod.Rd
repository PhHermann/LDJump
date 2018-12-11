\name{calcRegMod}
\alias{calcRegMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Regression Model under User Input Demography (Scenario)
%%  ~~function to do ... ~~
}
\description{This function computes the regression model for user input demographic scenarios. Moreover, the user is able to handle the sample sizes, lengths, and recombination rates of the simulated populations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcRegMod(n = c(10,16,20), len = c(500,1000,2000,3000,5000), thth = 0.01, nsim = 100,
           fr = c(), pathToScrm, scenario, pathToMs2dna, status = T, pathLDhat, pathPhi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{A numeric vector containig by default 10, 16, and 20 reflecting the sample sizes of the simulated populations. It can be adapted to any vector. }
  \item{len}{A numeric vector containing the lengths of simulated sequences of the populations. By default 0.5, 1, 2, 3, and 5 kb but can be adapted to any integer values.}
  \item{thth}{A numeric value for the mutation rate theta under which the populations are simulated. By default 0.01 but can be adapted to any numeric value.}
  \item{fr}{A numeric vector containing the recombination rates under which one wants to simulate. By default it is set to an empty vector and uniform random variables are simulated from 5 intervals with \code{nsim} values per interval.}
  \item{nsim}{An integer value for the number of replications (populations) simulated per setup. Setups result from all combinations of sample sizes and sequence lengths. This value can be adapted to any integer value.}
  \item{pathToScrm}{
  A character string containing the path to scrm. This path and the installation of \href{https://github.com/scrm/scrm/wiki/Command-Line-Options}{scrm} is necessary for the computation of the function.  }
  \item{scenario}{A character string containing the demography model (scenario) under which the populations should be simulated. We refer to  \href{https://github.com/scrm/scrm/wiki/Command-Line-Options}{scrm} for details on how to define varying population sizes using the simulation package \code{scrm}.}
  \item{pathToMs2dna}{
  A character string containing the path to \code{ms2dna}. This path and the installation of \href{http://guanine.evolbio.mpg.de/bioBox/ms2dna_1.16.tgz}{ms2dna} is necessary for the computation of the function.}
  \item{status}{
  an optional logical value: by default \code{TRUE} such that the current processing status of the number of simulated populations is printed.
  }
    \item{pathLDhat}{
  A character string containing the path to LDhat. This path and the installation of \href{https://github.com/auton1/LDhat}{LDhat} is necessary for the computation of the package.
  }
  \item{pathPhi}{
  A character string containing the path to PhiPack. This path and the installation of \href{https://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz}{PhiPack} is necessary for the computation of the package.
  }
%%     ~~Describe \code{x} here~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{regMod}{The generalized additive regression model on the box-cox transformed true recombination rates using computed summary statistics from simulated populations under a user defined demography (scenario).}
  \item{data.all}{A data-frame containing all summary statistics per column and simulated samples of populations per row.}
%% ...
}
\references{
Auton, A. and McVean, G. (2007). Recombination rate estimation in the presence
of hotspots. Genome Research, 17(8), 1219–1227.

Frick, K., Munk, A., and Sieling, H. (2014). Multiscale change-point inference.
Journal of the Royal Statistical Society: Series B, 76(3), 495–580.

Futschik, A., Hotz, T., Munk, A., and Sieling, H. (2014). Multiscale DNA
partitioning: Statistical evidence for segments. Bioinformatics, 30(16), 2255–2262.

Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new   tools for the analysis of
genome-wide SNP data. Bioinformatics. doi: 10.1093/bioinformatics/btr521

Knaus BJ and Grünwald NJ (2017). VCFR: a package to manipulate and visualize
variant call format data in R. Molecular Ecology Resources, 17(1), pp. 44-53.
ISSN 757, <URL:http://dx.doi.org/10.1111/1755-0998.12549>.

McVean, G. A. T., Myers, S. R., Hunt, S., Deloukas, P., Bentley, D. R., and Donnelly,
P. (2004). The fine-scale structure of recombination rate variation in the human
genome. Science, 304(5670), 581–584.

Paradis E., Claude J. & Strimmer K. 2004. APE:  analyses of phylogenetics and evolution
in R language. Bioinformatics 20: 289-290.

The 1000 Genomes Project Consortium (2015). Aglobal reference for human genetic
variation. Nature, 526(7571), 68–74.

Wood, S.N. (2011) Fast stable restricted maximum   likelihood and marginal
likelihood estimation of semiparametric generalized linear models.
Journal of the Royal Statistical Society (B) 73(1):3-36

}
\author{
Philipp Hermann \email{philipp.hermann@jku.at}, Andreas Futschik
}

\note{
  This function only works with unix and having \href{http://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz}{PhiPack} installed. Optionally when also having \href{https://github.com/auton1/LDhat}{LDhat} (Auton and McVean (2007)) installed \code{LDJump} will compute estimates much faster. Hence, please properly check all paths to \href{http://www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz}{PhiPack} and in case also \href{https://github.com/auton1/LDhat}{LDhat} as well as the sequence files.
Moreover, the software packages \href{https://github.com/scrm/scrm/wiki/Command-Line-Options}{scrm} and \href{http://guanine.evolbio.mpg.de/bioBox/ms2dna_1.16.tgz}{ms2dna} need to be installed for simulating populations under a user input demography (scenario).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{link{LDJump}}, \code{\link{summary_statistics}}, \code{\link{vcfR_to_fasta}}, \code{\link{getPhi}}, \code{\link{get_smuce}}, \code{\link[stepR]{smuceR}}, \code{\link[quantreg]{rq}}, \code{\link[mgcv]{gam}}, \code{\link[vcfR]{vcfR2DNAbin}}, \code{\link[genetics]{diseq}}, \code{\link[genetics]{genotype}}, \code{\link[Biostrings]{readDNAStringSet}}
}

\examples{
##### Do not run these examples                                         #####
##### scenario =  " -eG 0.0 0 -eG 0.42 -100 -eG 0.5 100 "               #####
##### simulatedData = calcRegMod(nsim=100,pathToScrm="/path/To/Scrm/",  #####
#####                scenario=scenario,pathToMs2dna="/path/To/Ms2dna/", #####
#####                 pathLDhat = "/path/to/LDhat/",                    #####
#####                 pathPhi = "/path/to/Phi/")                        #####
##### regMod = simulatedData[[1]]                                       #####
##### result = LDJump(fileName, alpha = 0.05, segLength = 1000,         #####
#####                 pathLDhat = "/path/to/LDhat/",                    #####
#####                 pathPhi = "/path/to/Phi/",                        #####
#####                 format = "fasta", regMod = regMod)                #####
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{models, regression}
\keyword{datasets}
\keyword{list}
\keyword{htest}
\keyword{datagen}
