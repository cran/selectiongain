\name{multistagetp}
\alias{multistagetp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the truncation point}
\description{
  
  This function calculates the coordinates of the truncation point Q for given selected fractions \code{alpha}. The R function \code{\link{uniroot}} in core package stats is called internally to solve the truncation point equations. 
}
\usage{
multistagetp(alpha, corx, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{Probability vector \eqn{\vec{\alpha}}, with length n, for random variable X. In other disciplines, such as plant breeding, it is also called the selected fraction.}  
  
  
  \item{corx}{n-dimensional matrix. It is the correlation matrix of selection indices X. More detail see \code{\link{multistagegain}}.}
    
\item{alg}{An object used to switch between two algorithms. More detail see \code{\link{multistagegain}}.}

}

\details{ 

This function calculates the non-equi coordinate quantile vector \eqn{Q=\{q_{1},q_{2},...,q_{n}\}} for a multivariate normal distribution from a given selection fraction vector \eqn{\textbf{A}=\{ \alpha_{1},\alpha_{2},...,\alpha_{n} \}}. It can be compared with the function \code{\link{qmvnorm}} in package \pkg{mvtnorm}, which calculates the equi-coordinate quantile  \eqn{q} for multi-variate normal distribution from a given alpha scalar.
 
 
The truncation point \eqn{\textbf{Q} } must be calculated from the vector of selected fractions \eqn{\textbf{A} }. The coordinates \eqn{q_i } are obtained by inverting the following equations sequentially: 

\deqn{ 
\alpha_1= f (q_1) =\int_{q_1}^{\infty} \phi_1 (x_1) \, dx_1 \label{equ-alpha1},
} 
which yields \eqn{q_1 };
 
\deqn{
\alpha_1 \alpha_2= f (q_1,q_2) =\int_{q_2}^{\infty} \int_{q_1}^{\infty}  \phi_2 (x_1,x_2; \rho_{x_1,x_2})\, dx_1 \,dx_2,
}

which yields \eqn{q_2 } for given values of \eqn{q_1 }, \eqn{\alpha_1 } and \eqn{\alpha_2 }, and so on for \eqn{i=3,4,\ldots,n }. Here we have assumed that \eqn{\bm{\Sigma} } is diagonal, so the above inverting function can be simplified into

\deqn{
 \alpha_i=\int_{q_i}^{\infty} \phi_1 (x_i)\, d x_i \label{equ-simple}
}

If \eqn{\bm{\Sigma} } is not diagonal, the value of \eqn{q_n } is finally obtained by using \eqn{q_{1},q_{2},\ldots,q_{n-1} } and \eqn{\alpha_{1},\alpha_{2},\ldots,\alpha_{n} }:

\deqn{
 \prod_1^n \alpha_i=\int_{q_n}^{\infty}\int_{q_{n-1}}^{\infty}\ldots\int_{q_1}^{\infty} \phi_n (\bm{X};\bm{\Sigma})\, d \textbf{x} \label{equ-alpha2}
}

These equations can be solved numerically via a fast root search algorithm, which was described for example by Brent (1973) and implemented in the function \code{\link{uniroot}}, though the function \code{\link{qmvnorm}}. 


}
\value{
  The out come is a vector for the coordinates.
}
\references{
Brent, R. (1973) Algorithms for Minimization without Derivatives. Prentice-Hall, Englewood Cliffs, New Jersey.
 }
\author{ Xuefei Mi}
\note{ When a vector \eqn{A} is given, the quantiles are calculated consecutively to satisfy the given \eqn{A}. The calculation from other direction to \eqn{-\infty} of the integral is also possible for \code{\link{qmvnorm}}. But not for this function, the user have to reverse the integral.  
}
\seealso{ \code{\link{multistagegain}}}
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
