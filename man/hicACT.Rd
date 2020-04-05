\name{hicACT}
\alias{hicACT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HiC-ACT: Aggregated Cauchy Test approach for improved chromatin interaction calling from Hi-C data
%%  ~~function to do ... ~~
}
\description{hicACT() computes HiC-ACT smoothed p-values of chromatin interaction calls from Hi-C data
}
\usage{
hicACT(infile, kb, h, pthres, outdir, outname="ACT_adjusted",
                   bin1col=2, bin2col=4, ccountcol=5, pcol=6,
                   ignore_warnings=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{output file from a Hi-C chromatin interaction calling method, such as Fit-Hi-C (column names required)}
  \item{kb}{data resolution in Kb, e.g. 10}
  \item{h}{smoothing parameter, see 'Details' for suggestions}
  \item{pthres}{p-value threshold for selecting which p-values to adjust using HiC-ACT}
  \item{outdir}{output directory path (defaults to working directory)}
  \item{outname}{desired output file name}
  \item{bin1col}{column number containing the first bin IDs (defaults to 2 for Fit-H-C)}
  \item{bin2col}{column number containing the second bin IDs (defaults to 4 for Fit-Hi-C)}
  \item{ccountcol}{column number containing the observed contact counts (defaults to 5 for Fit-Hi-C)}
    \item{pcol}{column number containing the p-values (defaults to 6 for Fit-Hi-C))}
    \item{ignore_warnings}{default is FALSE}
}
\details{
Suggestions for the smoothing parameter (h) based on data resolution:

  - If 5Kb, then h=40.

  - If 10Kb, then h=20.

  - If 20Kb, then h=12.

  - If 25Kb, then h=11.

  - If 40Kb, then h=5.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
