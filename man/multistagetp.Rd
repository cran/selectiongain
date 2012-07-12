\name{multistagetp}
\alias{multistagetp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the truncation point}
\description{
  
  This function calculates the coordinates of the truncation point Q for given selected fractions A. The R function uniroot in core package stats is called internally to solve the truncation point equations. 
}
\usage{
multistagetp(alpha, corx, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{Probability vector \eqn{\vec{\alpha}} for random variable X. In plant breeding, it is also called the selected fraction.}  
  \item{corx}{Correlation matrix of X. }
  \item{alg}{It decides which algorithm will be used, the Genz and Bretz' algorithm is used by default. The Miwa's algorithm can be optional (Mi. et al. 2009; Genz. et al. 2010).}
}
\details{ 

This function calculates the non-equi coordinate quantile vector \eqn{Q=\{q_{1},q_{2},...,q_{n}\}} for a multivariate normal distribution from a given alpha vector \eqn{\vec{\alpha}=\{ \alpha_{1},\alpha_{2},...,\alpha_{n} \}}. It can be compared with the function qmvnorm() in R-package \pkg{mvtnorm}, which calculates the equi coordinate quantile  \eqn{q} for multi-variate normal distribution from a given alpha scalar.
 
More details are in the JSS paper section 3.2.

}
\value{
  The out come is a vector for the coordinates.
}
\references{

A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9, 2010.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470, 1951

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

X. Mi, H.F. Utz. and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)
 }
\author{ Xuefei Mi}
\note{ When a vector \eqn{A} is given, the quantiles are calculated consecutively to satisfy the given \eqn{A}. The calculation from other direction to \eqn{-\infty} of the integral is also possible for qmvnorm(). But not for this function, the user have to reverse the integral.  
}
\seealso{ selectiongain(), qnorm() }
\examples{
 dim=3
 alpha=rep(0,dim)

 alpha[1] = 0.3333
 alpha[2] = 0.25
 alpha[3] = 0.2

corr=diag(3)

corr12=0.3016
corr[1,2]=corr12
corr[2,1]=corr12

corr13=0.5630
corr[1,3]=corr13
corr[3,1]=corr13

corr23=0.5630
corr[2,3]=corr23
corr[3,2]=corr23

multistagetp(alpha,corx=corr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal}
% __ONLY ONE__ keyword per line
