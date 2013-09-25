\name{multistagetp}
\alias{multistagetp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the truncation points}
\description{ This function calculates the coordinates of the truncation points Q for given selected fractions \eqn{\vec{\alpha}=\{ \alpha_{1},\alpha_{2},...,\alpha_{n} \}} and correlation matrix of X. The R function uniroot in core package stats is called internally to solve the truncation point equations. It is needed for multistagegain function to calculated the expected selection gain.
}
\usage{
multistagetp(alpha,  corr,  alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{(n-dimensional vector) is probability vector \eqn{\vec{\alpha}} for random variable X. In plant breeding, it is also called the selected fraction.}  
     
  \item{corr}{(n+1-dimensional matrix) is the correlation matrix of y and X, which is introduced in function multistagecorr. The correlation matrix must be symmetric and positive-definite. The correlation matrix estimated from practices sometimes can be negative-definite, it must be adjusted before using this function. Before starting the calculations, the user is recommended to check the correlation matrix. The default value (in this function only) is \code{NA}. }


  
\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz et al. (2009, 2013), will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, its computational speed is slower. We recommend to use Miwa algorithm of this parameter.}

 }

\details{ 

This function calculates the non-equi coordinate quantile vector \eqn{Q=\{q_{1},q_{2},...,q_{n}\}} for a multivariate normal distribution from a given \eqn{\vec{\alpha}}. It can be compared with the function qmvnorm() in R-package \pkg{mvtnorm}, which calculates only the equi coordinate quantile  \eqn{q} for multi-variate normal distribution from a given alpha scalar.
 
}
\value{
  The output is a vector of the coordinates.
}
\references{

A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9995, 2013.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

 }
\author{ Xuefei Mi}
\note{ When a vector \eqn{A} is given, the quantiles are calculated consecutively to satisfy the given \eqn{A}. The calculation from other direction to \eqn{-\infty} of the integral is also possible for qmvnorm(). But not for this function, the user have to reverse the integral.  
}
\seealso{ selectiongain(), qnorm() }
\examples{

# first example

VCGCAandError=c(0.40,0.20,0.20,0.40,2.00)
VCSCA=c(0.20,0.10,0.10,0.20)

corr.matrix = multistagecor(maseff=0.40, VGCAandE=VCGCAandError,
VSCA=VCSCA, T=c(1,1,5), L=c(1,3,8), Rep=c(1,1,1))

N1=4500;N2=919;N3=45;Nf=10

Q=multistagetp(c(N2/N1,N3/N2,Nf/N3),  corr=corr.matrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal}
% __ONLY ONE__ keyword per line
