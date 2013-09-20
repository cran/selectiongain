 \name{multistageoptimum.search}
\alias{multistageoptimum.search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing three-stage selection in plant breeding with one marker-assisted selection stage and two phenotypic selection stages}

\description{ This function is used to calculate the maximum of \eqn{\Delta G} with correlation matrix, which is depending on locations, testers and replicates by grid search algorithm. The changing correlation matrix of three-stage selection are the testcross progenies of DH lines in one marker-assisted selection (MAS) stage and two phenotypic selection (PS) stages.
}
\usage{
multistageoptimum.search (maseff=0.4, VGCAandE, 
  VSCA, CostProd, CostTest,  Nf, Budget, N2grid, 
  N3grid, L2grid, L3grid, T2grid, T3grid, R2, R3, alg, detail, fig)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{maseff}{ is the efficiency of MAS.}
   
\item{VGCAandE}{is the vector of variance components for general combining ability and the plot error.}  
 
 \item{VSCA}{is the vector of variance components for specific combining ability.}  

\item{CostProd}{ contains the initial costs of producing or identifying a candidate in each stage.}

\item{CostTest}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}

\item{Nf}{ is the number of final selected candidates.}


\item{Budget}{contains the value of total budget.}


 \item{N2grid}{is the vector of lower and upper limits as well as the grid width of number of candidates in the first field test stage.}

\item{N3grid}{is the vector of lower and upper limits as well as the grid width of number of candidates in the second field test stage.}

\item{L2grid}{is the vector of lower and upper limits of number of location as well as the width in the first field test stage.}

\item{L3grid}{is the vector of lower and upper limits of number of location as well as the width in the second field test stage.}

\item{T2grid}{is the vector of lower and upper limits of number of tester as well as the width in the first field test stage.}

\item{T3grid}{is the vector of lower and upper limits of number of tester as well as the width in the second field test stage.}

\item{R2}{is the number of replicates in the first field test stage. By default it is 1. }

\item{R3}{is the number of replicates in the second field test stage. By default it is 1.}

\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz et al. (2009, 2013), will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, its computational speed is slower. We recommend to use Miwa algorithm of this parameter.}

\item{detail}{ is the control parameter to decide if the result of all the grids will be given (\code{=TRUE}) or only the maximum (\code{=FALSE}).}


\item{fig}{is the control parameter to decide if a figure of contour plot will be saved in the default folder. The default value is \code{FALSE}, which means no figure will be saved.}

}
%\details{ 

%More details are in the JSS paper section 3.5.
 
%}
\value{
If \eqn{\texttt{detail}} = FALSE, the output of this function is a vector of the optimum allocation i.e., which achieves the maximum \eqn{\Delta G}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table.}
\references{

A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9995, 2013.

E.L. Heffner, A.J. Lorenz, J.L. Jannink, and M.E. Sorrells.  Plant breeding with genomic selection: gain per unit time and cost. Crop Sci. 2007 50: 1681-1690.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.



 }
\author{ Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{

CostProd =c(0.5,1,1)
CostTest = c(0.5,1,1)
Budget=1021
# Budget is very small here to save time in package checking
# for the example in Heffner's paper, please change it to Budget=10021

VCGCAandError=c(0.4,0.2,0.2,0.4,2)
VCSCA=c(0.2,0.1,0.1,0.2)
Nf=10


multistageoptimum.search (maseff=0.4, VGCAandE=VCGCAandError, 
VSCA=VCSCA, CostProd = c(0.5,1,1), CostTest = c(0.5,1,1), 
Nf = 10, Budget = Budget, N2grid = c(11, 1211, 30), 
N3grid = c(11, 211, 5), L2grid=c(1,1,1), L3grid=c(6,6,1),
T2grid=c(1,2,1), T3grid=c(3,5,1), R2=1, R3=1, alg = Miwa(), 
detail=TRUE, fig=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
