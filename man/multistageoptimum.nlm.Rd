\name{multistageoptimum.nlm}
\alias{multistageoptimum.nlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing n-stage selection with the NLM algorithm for a given correlation matrix}
\description{This function is used to calculate the maximum of \eqn{\Delta G} with given correlation matrix by non-linear minimization algorithm.
  }
\usage{
multistageoptimum.nlm(corr, ini.value, 
Budget, CostC, CostTv, Nf, iterlim, alg,N.upper, N.lower)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{corr}{is the correlation matrix of y and X, which is introduced in function multistagecorr. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix.}

\item{ini.value}{is a vector, which stores the number of candidates in each stage for the algorithm to begin with. As default, it will use \eqn{N=\{N_1,N_2,...,N_n\}=\{a+1,...,a+n\}}, where a is defined as \eqn{Budget/(CostC+sum(CostTv)+1)}}.


\item{Budget}{contains the value of total budget.}

\item{CostC}{ contains the costs of producing or identifying a candidate.}

\item{CostTv}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}

\item{Nf}{ is the number of final selected candidates.}

\item{iterlim}{is the maximum number of iterations to be executed before the Newton algorithm is terminated. By default it is equal to 20. If the \eqn{\texttt{Budget}} increases 10 times for making the selection, the value of \eqn{\texttt{iterlim}} has to be increased \eqn{lg(10)} times.}

\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz (2009) will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et. al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, the computational speed is slower. We recommend the user to use Miwa algorithm of this parameter.}

\item{N.upper}{ is the vector of up limits of number of candidates X.}

\item{N.lower}{is the vector of low limits of number of candidates X.}



}
%\details{

%More details are in the JSS paper section 3.6.
 
%}
\value{
The output of this function is a vector similar as in multistageoptimal.grid() if detail = FALSE. However, the optimal number of candidates in each stage determined by the NLM algorithm is clearly not an integer, because the function uses a numerical algorithm, which depends on derivatives. 
}
\references{
A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9, 2010.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

X. Mi, H.F. Utz. and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470, 1951.

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

multistageoptimum.nlm(N.upper=rep(100,3), corr=corr, Budget=200,
 CostC=0.5, Nf=5)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
