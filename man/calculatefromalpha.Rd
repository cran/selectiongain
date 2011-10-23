\name{calculatefromalpha}
\alias{calculatefromalpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Non-equi coordinate quantile for singly-truncated multivariate normal distribution}
\description{This function calculates the non-equi coordinate quantile vector \eqn{Q=\{q_{1},q_{2},...,q_{n}\}} for a multivariate normal distribution from a given alpha vector \eqn{\vec{\alpha}=\{ \alpha_{1},\alpha_{2},...,\alpha_{n} \}}. It can be compared with the function qmvnorm() in R-package \pkg{mvtnorm}, which calculates the equi coordinate quantile  \eqn{q} for multi-variate normal distribution from a given alpha scalar.

E.g. Let \eqn{f} be the density function, for  \eqn{n=2}:

  the equi coordinate quantile \eqn{\vec{q}=\{q,q,...,q\} } defined by \pkg{mvtnorm} is:

  \eqn{ \alpha=\int_{q}^{\infty}\int_{q}^{\infty} f(x_1,x_2) dx_1 dx_2}
  
  the non-equi coordinate quantile \eqn{Q=\{q_{1},q_{2},...,q_{n}\} } defined by us is:
  
  \eqn{\alpha_1=\int_{q_1}^{\infty} f(x_1) dx_1}
  
  \eqn{ \alpha_2=\int_{q_1}^{\infty}\int_{q_2}^{\infty} f(x_1,x_2) dx_1 dx_2}
 
  
  
}
\usage{
calculatefromalpha(alpha, dim, corr, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{Dimension of the random variable X.}
  \item{alpha}{Probability vector \eqn{\vec{\alpha}} for random variable X. In plant breeding, it is also called the selected fraction.}  
 
  \item{corr}{ Correlation matrix }
  \item{alg}{It decides which algorithm will be used, the Genz and Bretz' algorithm is used by default. The Miwa's algorithm can be optional (Mi. et al. 2009; Genz. et al. 2010)}


}
\details{ There is a limitation, if the absolute value of \eqn{q_i} is greater than 100, then the absolute value is settled to be 100.
 

}
\value{
  The out come is a vector for the coordinates.
}
\references{

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

A., Genz, F., Bretz (2009). Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg.

A., Genz, F., Bretz, T., Miwa, X., Mi, F., Leisch, F., Scheipl, T., Hothorn, 2010. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9.



 }
\author{ Xuefei Mi}
\note{ When a vector \eqn{\vec{\alpha}} is given, the quantiles are calculated consecutively to satisfy the given \eqn{\vec{\alpha}}.  The calculation from other direction to minus infinite of the integral is also possible for qmvnorm(). But not for this function, the user have to reverse the integral.  

}
\seealso{ selectiongain(), qnorm() }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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

calculatefromalpha(alpha,dim=3,corr)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Non-equi coordinate quantile (NCQ) }
% __ONLY ONE__ keyword per line
