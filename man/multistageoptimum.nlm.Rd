\name{multistageoptimum.nlm}
\alias{multistageoptimum.nlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing n-stage selection with the NLM algorithm for a given correlation matrix}
\description{This function is used to calculate the maximum of \eqn{\Delta G} with given correlation matrix by non-linear minimization algorithm.
  }
\usage{
multistageoptimum.nlm(corr, Vg, ini.value, 
Budget, CostProd, CostTest, 
Nf, iterlim, alg, N.upper, N.lower)

}



%- maybe also 'usage' for other objects documented here.
\arguments{

\item{corr}{is the correlation matrix of y and X, which is introduced in function multistagecorr. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix.}

 \item{Vg}{is genotypic variance \eqn{\delta_y^2}. The default value is 1.}

\item{ini.value}{is a vector, which stores the number of candidates in each stage for the algorithm to begin with. As default, it will use \eqn{N=\{N_1,N_2,...,N_n\}=\{a+1,...,a+n\}}, where a is defined as \eqn{(N.upper+N.lower)/4}}.

\item{Budget}{contains the value of total budget.}

\item{CostProd}{ contains the initial costs of producing or providing a candidate in each stage}

\item{CostTest}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}


\item{Nf}{ is the number of finally selected candidates.}

\item{iterlim}{is the maximum number of iterations to be executed before the Newton algorithm is terminated. By default it is equal to 20. If the \eqn{\texttt{Budget}} increases 10 times for making the selection, the value of \eqn{\texttt{iterlim}} has to be increased \eqn{lg(10)} times.}

\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz et al. (2009, 2013), will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, its computational speed is slower. We recommend the user to use Miwa algorithm of this parameter.}

\item{N.upper}{ is the vector of up limits of number of candidates X.}

\item{N.lower}{is the vector of low limits of number of candidates X.}



}
%\details{

%More details are in the JSS paper section 3.6.
 
%}
\value{
The output of this function is a vector similar as in multistageoptimal.grid(). However, the optimal number of candidates in each stage determined by the NLM algorithm is clearly not an integer, because the function uses a numerical algorithm, which depends on derivatives. 
}
\references{
A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9995, 2013.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. J. Royal Stat. Soc., Ser. B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselektion in der Pflanzenzuechtung (in German). Doctor thesis, University Hohenheim, 1969.

W.G. Cochran. Improvement by means of selection. In J. Neyman (ed.) Proc. 2nd Berkeley Symp. on Mathematical Statistics and Probability. University of California Press, Berkeley., 1951.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

 }
\author{ Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{

 VCGCAandError=c(0.40,0.20,0.20,0.40,2.00)
 VCSCA=c(0.20,0.10,0.10,0.20)

corr = multistagecor (maseff=0.40,
  VGCAandE=VCGCAandError,  VSCA=VCSCA, T=c(1,1,5),
  L=c(1,3,8), Rep=c(1,1,1))

# the time of nlm have to be controled in 5 s, so this example will not be uploaded into cran

#multistageoptimum.nlm( corr=corr, Vg=0.4,
#Budget=1021, CostProd=c(0.5,0,0),CostTest=c(0.5,6,40), Nf=10,
# N.upper=c(600,120,20), N.lower=rep(5,3))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
