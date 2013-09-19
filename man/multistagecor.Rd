\name{multistagecor}
\alias{multistagecor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating correlation matrix in a plant breeding context}
\description{
This function is used to calculate the correlation matrix \eqn{\bm{\Sigma}^{*}} of y and X (where y is the true genotypic value and \eqn{\bm{X}=\{X_1,...X_n \}} are the selection indices from each selection stage) for n-stage selection in plant breeding. It is assumed that the genetic structure of the candidates to be selected are genetically fixed, e.g., potential cultivars, clones, inbred lines or testcross progenies of inbred lines with the same or different testers in all stages.

  
}
\usage{
multistagecor(maseff,VGCAandE,VSCA,T,L,Rep,index,covscatype,detail)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{maseff}{ is the efficiency of Marker-Assisted Selection (MAS). The default value is NA. If a value between 0 and 1 is assigned to \code{maseff}, then the first selection stage will be considered as MAS (Heffner 2010).}
   
  \item{VGCAandE}{is the vector of variance components of genetic effect, gene by location interaction, gene by year interaction, gene by location  by year interaction and the plot error. When \code{VSCA} is specified, it stands for the general combining ability.  The default value is {1,1,1,1,1}. Variances types listed in Longin (2007) can be used. E.g., \code{VGCAandE="VC2"} will set the value as {1,0.5,0.5,1,2}. }  
 
  \item{VSCA}{is the vector of variance components for specific combining ability. The default value is {0,0,0,0}. }  
 
 \item{T}{is the vector of number of testers at each stage. If there is no testers applied in certain stage, the value of testers at this stage have to be 1. }  

  \item{L}{is the vector of number of locations at each stages.}
  
  \item{Rep}{is the vector of number of replicates at each stage.}  
 
  \item{index}{is the control parameter. If it equals TRUE, the optimum selection index of the Longin (2007) will be used in the calculation of correlation matrix without MAS. 
 }

   \item{covscatype}{is the type of the covariance. Longin's type (\code{covscatype}=c("LonginII")) is by default. }
   

   \item{detail}{ is the control parameter to decide if the correlation matrix, optimal selection index and covariance matrix will be given or only the correlation matrix. The default value of it is FALSE.}

}
\value{
The default output is a matrix with dimension n+1 and can be used as input parameter of function multistagegain. When  value of detail=TRUE, the result of correlation matrix, optimal selection index and covariance matrix will be given.
}
\references{

C. Longin, H.F. Utz., J. Reif, T. Wegenast, W. Schipprack and A.E. Melchinger. Hybrid maize breeding with doubled haploids: III. Efficiency of early testing prior to doubled haploid production in two-stage selection for testcross performance. Theor. Appl. Genet. 2007 115: 519-527.

E.L. Heffner, A.J. Lorenz, J.L. Jannink, and M.E. Sorrells.  Plant breeding with genomic selection: gain per unit time and cost. Crop Sci. 2007 50: 1681-1690.


 }
\author{ Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{
# example for calculating correlation matrix without MAS 
multistagecor(VGCAandE=c(1,0.5,0.5,1,2),L=c(2,10),T=c(1,1),Rep=c(1,1))
multistagecor(VGCAandE="VC2",L=c(2,10),T=c(1,1),Rep=c(1,1),index=TRUE)

# example for calculating correlation matrix with MAS in the first stage 
VCgca=c(0.40,0.20,0.20,0.40,2.00)
VCsca=c(0.20,0.10,0.10,0.20)
corr.matrix = multistagecor (maseff=0.40, VGCAandE=VCgca,
  VSCA=VCsca, T=c(1,1,5), L=c(1,3,8), Rep=c(1,1,1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
