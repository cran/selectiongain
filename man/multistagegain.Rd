\name{multistagegain}
\alias{multistagegain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the multi-stage selection gain}
\description{
 
This is the main function of the package and uses following equation given by Tallis(1965):

\eqn{
\frac{\partial m(\textbf{t})}{\partial t_i}|_{\textbf{t}=\textbf{0}}=  E(X_i) =\frac{1}{\alpha} \sum_{k=0}^{n} \rho_{i,k} \,  \phi_1(q_k) \,  \Phi_{n} (A_{k,s};R_k) 
}

 to calculate selection gain for given correlation matrix and coordinates of the truncation points. 



}
\usage{
multistagegain(Q, corr, alg, lim.y, stages=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 
\item{Q}{ (length n) refers to the coordinates of the truncation points Q, which is the output of the next function (\pkg{multistagetp}) that we are going to introduced.}
  
\item{corr}{(n+1-dimensional matrix) is the correlation matrix of y and X. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix, which is usually obtained by analysis of one or several experiments designed for estimating covariance components as well as phenotypic and genotypic correlations.}

\item{alg}{is used to switch between two algorithms. If alg = GenzBretz(), which is by default, the quasi-Monte Carlo algorithm from Genz(1999) will be used. If alg = Miwa(), the program will use the Miwa algorithm (Mi2009), which an analytical solution of the MVN integral (Miwa2003). Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits), however the computational speed is also slower. Therefore, we recommend the user to use the default value of this parameter.}
 
\item{lim.y}{ is the lower limit of y and is set to -200 as default, which is on the safe side.} 
  
\item{stages}{If = TRUE, the gain from each stages will be shown. }



}
\details{

This function calculates the well-known selection gain \eqn{\Delta G(y)}, which is described by Cochran (1951). For one-stage selection the gain is defined as \eqn{\Delta G(y) = i \rho_{y} \rho_{1}}, where \eqn{i} is the selection intensity, \eqn{\rho_{1}} is the correlation between the true breeding value and the selection index \eqn{y} (Utz 1969). 

More details are in the JSS paper section 3.1.



}

\value{
  The value returned, is the gain of selection.
  }
\references{


A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9, 2010.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470, 1951.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

X. Mi, H.F. Utz. and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)

}


\author{ Xuefei Mi
}
\note{ No further notes

}
\seealso{
 No link
}
\examples{

k=c(-200,0.4308,0.9804,1.8603)
Q=c(0.4308,0.9804,1.8603)
corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)



multistagegain(Q=Q,corr=corr)


#####
# code for testing the parameter stages
#####

multistagegain(Q=Q,corr=corr,stages=TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal }

