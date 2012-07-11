\name{multistageoptimal.grid}
\alias{multistageoptimal.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing multi-stage selection with grid algorithm for fixed correlation matrix}
\description{This function is used to calculate the maximum of \eqn{\Delta G(y)} with fixed correlation matrix.
  
}
\usage{
multistageoptimal.grid(N.upper, N.lower=rep(1,length(N.upper)), corr, 
num.grid=11, Budget, CostC=1, CostTv=rep(1,length(N.upper)),N.fs, 
detail=FALSE,alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 
\item{N.upper}{is the vector of up limits of number of candidates X.}

\item{N.lower}{is the vector of low limits of number of candidates X.}


\item{corr}{Correlation matrix of X. }

\item{num.grid}{ is the number of equally distanced points which divided the axis of \eqn{x_1} into \eqn{num.grid-1} intervals and there are \eqn{(num.grid-1)^n} grids in a n dimensional hyper cube. The default value of num.grid is 11, so the hyper cube is divided into \eqn{10^n} parts by the grid points. If the grid point is more than the maximum number of certain stage i, then the grid number of this stage will equal to the number of grid points, i.e., if \eqn{num.grid > N_i}, then the number of grid points for the i-th axis is \eqn{N_i}.}

\item{Budget}{contains the value of total budget.}

\item{CostC}{ contains the costs of producing or identifying a candidate.}

\item{CostTv}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}

\item{N.fs}{ is the number of final selected candidates.}

\item{detail}{ is the control parameter to decide if the result of all the grids will be given or only the maximum.}

\item{alg}{is used to switch between two algorithms. If alg = GenzBretz(), which is by default, the quasi-Monte Carlo algorithm from Genz(1999) will be used. If alg = Miwa(), the program will use the Miwa algorithm (Mi2009), which an analytical solution of the MVN integral (Miwa2003). Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits), however the computational speed is also slower. Therefore, we recommend the user to use the default value of this parameter.}

}
\details{ 

More details are in the JSS paper section 3.5.
 

}
\value{
If \eqn{\texttt{detail}} = FALSE, the output of this functions is a vector with the optimal number of candidates in each stage (\eqn{\textbf{N}}) and the maximum \eqn{\Delta G(y)}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table.}
\references{

A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9, 2010.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470, 1951.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

X. Mi, H.F. Utz. and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)


 }
\author{ Xuefei Mi}
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

multistageoptimal.grid(N.upper=rep(100,3), corr=corr, Budget=200, CostC=0.5, N.fs=5)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
