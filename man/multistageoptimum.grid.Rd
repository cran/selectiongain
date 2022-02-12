\name{multistageoptimum.grid}
\alias{multistageoptimum.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing multi-stage selection with grid algorithm for a given correlation matrix}
\description{This function is used to calculate the maximum of \eqn{\Delta G} for a given correlation matrix by grid search algorithm.
  
}
\usage{
multistageoptimum.grid(corr, Vg,
num.grid, width, Budget, CostProd,
CostTest,Nf,alg,detail,fig,N.upper, N.lower,alpha.nursery,cost.nursery,vargain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 \item{Vg}{is genotypic variance \eqn{\delta_y^2}. The default value is 1.}
 
\item{corr}{is the correlation matrix of y and X, which is introduced in the function multistagecorr. The correlation matrix must be symmetric and positive-definite. If the estimated correlation matrix is negative-definite, it must be adjusted before using this function. Before starting the calculations, it is recommended to check the correlation matrix.}


\item{num.grid}{ is the number of equally distanced points that divided the axis of \eqn{x_1} into \eqn{num.grid_i-1} intervals and there are \eqn{\prod_i(num.grid_i)} grids in a n dimensional hyper cube.   If \eqn{num.grid > N_i}, then the number of grid points for the i-th axis is \eqn{N_i}. The default value of it is \code{NA}.}

\item{width}{ is the width between the equally distanced points. The default value is \code{NA}.}


\item{Budget}{contains the value of total budget.}

\item{CostProd}{ contains the initial costs of producing or providing a candidate in each stage}

\item{CostTest}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}

\item{Nf}{ is the number of finally selected candidates.}

\item{detail}{ is the control parameter to decide if the result of all the grids will be given or only the maximum. The default value is \code{FALSE}.}


\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz et al. (2009, 2013), will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, its computational speed is slower. We recommend to use Miwa algorithm of this parameter.}
 


\item{fig}{is the control parameter to decide if a figure of contour plot will be saved in the default folder of R. The default value is \code{FALSE}, which means no figure will be saved.}

\item{N.upper}{is the vector of upper limits of number of candidates X.}

\item{N.lower}{is the vector of lower limits of number of candidates X.}

\item{alpha.nursery}{a value that should be 0<x<1, prelimitery test alpha fraction should be used for the stage 1.  it is setted to 1 as default, when no prelimitery test "nursery stage". }

\item{cost.nursery}{a vector of length two c([cost of producing a DH line],[cost of testing a DH in nursery]). The default value is {0,0}. }

\item{vargain}{is the logical variable to calculate the variance after multi-stage selection. Default is FALSE. Please see more details in the documentation for the function multistagevariance.The default value is \code{FALSE}}


}
\details{ 


for the new added to parameters "alpha.nursery" and "cost.nursery" since v2.0.47: 
 
 After producing new DH lines, breeders do NOT go directly for a selection stage in the field, neither for genomic selection. Most of the times, they prefer to make a small field experiment (called "nursery") in which all DH lines are observed and discarded for other traits as disease resistance. That means, all DH lines with poor resistance will be discarded. At the end of the nursery stage only certain amount of DH lines (alpha) advance to the first selection stage (phenotypic or genomic). Specially in maize that makes sense, because in experience around 90 percent of the new DH lines are very weak in terms of per se performance what make them not suitable as new hybrid parents. Then, budget should not be used to make genotyping on or testcrossing with them. Only the alpha fraction should be used for entering the stage 1 of the multistageoptimum.search function.

More details are available in the Crop Science and Computational Statistics papers.
 
 
}
\value{
If \eqn{\texttt{detail}} = FALSE, the output of this functions is a vector with the optimal number of candidates in each stage (\eqn{\textbf{N}}) and the maximum \eqn{\Delta G}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table.}
\references{


A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9995, 2013.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. J. Royal Stat. Soc., Ser. B, 23(1):223-229, 1961.

W.G. Cochran. Improvement by means of selection. In J. Neyman (ed.) Proc. 2nd Berkeley Symp. on Mathematical Statistics and Probability. University of California Press, Berkeley, 1951.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution. R Journal, 1:37-39, 2009.



 }
\author{ Xuefei Mi, Jose Marulanda}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{
corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)

Budget=200

multistageoptimum.grid( Vg=1, num.grid=11, corr=corr, Budget=Budget,
 CostProd=c(0.5,0,0), CostTest=c(0.5,1,1), Nf=5, 
N.upper=rep(Budget,3), N.lower=rep(1,3))

multistageoptimum.grid( Vg=1, num.grid=11, corr=corr, Budget=Budget,
 CostProd=c(0.5,0,0), CostTest=c(0.5,1,1), Nf=5, 
N.upper=rep(Budget,3), N.lower=rep(1,3),detail=TRUE,fig=TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
