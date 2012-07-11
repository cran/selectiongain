\name{multistagecor}
\alias{multistagecor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating correlation matrix in a plant breeding context}
\description{
This function is used to calculate the correlation matrix \eqn{\bm{\Sigma}^{*}} for n-stage selection in plant breeding. It is assumed that the genetic structure of the candidates to be selected are genetically fixed, e.g., potential cultivars, clones, inbred lines or testcross progenies of inbred lines with the same or different testers in all stages.

  
}
\usage{
multistagecor(V, L,Rep,index, coe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{V}{is the vector of variance components.}  
 
  \item{L}{is the vector of locations at each stages.}
  
  \item{Rep}{  is the vector of replicates R at each stage. The values are all equal to one by default.
 }  
 
  \item{index}{is the control parameter. If it equals TRUE, the optimal selection index (longin2007) will be used in the calculation. 
 }

  \item{coe}{is the control parameter. If it equals TRUE (as default) and index equals TRUE, the coefficient of the optimal selection index will be given as an output. In other situations, there will be no output. 
}

}
\details{ 
 More details are in the JSS paper section 3.7.

}
\value{
The output is a matrix with dimension n+1 and can be used as input parameter of function multistagegain.
}
\references{

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470, 1951

X. Mi, H.F. Utz and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)

C. Longin, H.F. Utz., J. Reif, T. Wegenast, W. Schipprack and A.E. Melchinger. Hybrid maize
breeding with doubled haploids: III. Efficiency of early testing prior to doubled haploid production in two-stage selection for testcross performance. Theoretical and Applied Genetics 2007 115(4): 519-527.



 }
\author{ Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{
multistagecor(V=c(1,0.5,0.5,1,2),L=c(2,10))
multistagecor(V="VC2",L=c(2,10),index=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
