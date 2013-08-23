\name{multistagegain}
\alias{multistagegain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the expected multi-stage selection gain}
\description{
 
This is the main function of the package and uses the following equation given by Tallis (1961) for y, which is the true genotypic value:

\eqn{
\frac{\partial m(\textbf{t})}{\partial t_0}|_{\textbf{t}=\textbf{0}}=  E(X_0=y) =\frac{1}{\alpha} \sum_{k=0}^{n} \rho_{0,k} \,  \phi_1(q_k) \,  \Phi_{n} (A_{k,s};R_k) 
}

 to calculate the expected selection gain defined by Cochran (1951) for given correlation matrix and coordinates of the truncation points. 

}
\usage{
multistagegain(corr, Q, Vg=1, alg, partial=FALSE)



}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{corr}{is the correlation matrix of y and X, which is introduced in function multistagecorr. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix.}

\item{Q}{are the coordinates of the truncation points, which are the output of the next function multistagetp that we are going to introduce.}
  
\item{Vg}{is true breeding value. The default value is 1.}

\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz (2009) will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et. al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, the computational speed is slower. We recommend the user to use Miwa algorithm of this parameter.}
 
  
\item{partial}{If \code{partial = TRUE}, the partial gain (\eqn{\rho_{0,k} \,  \phi_1(q_k) \,  \Phi_{n} (A_{k,s};R_k)}) from deleting stage k will be shown. The default value is \code{FALSE}.}

%\item{absepsmvtnorm}{is passed to the \code{absep} variable in \pkg{mvtnorm} for controlling the error. The default is 0.001, which is same as the value in \pkg{mvtnorm}.}


}
\details{

This function calculates the well-known selection gain \eqn{\Delta G}, which is described by Cochran (1951). For one-stage selection the gain is defined as \eqn{\Delta G = i V_g \rho_{1}}, where \eqn{i} is the selection intensity, \eqn{\rho_{1}} is the correlation between the true breeding value (its variance is \eqn{V_g}) and the selection index (Utz 1969). 

}

\value{
  The value returned, is the expected gain of selection.
  }
\references{
A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9995, 2013.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselektion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470, 1951.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

}


\author{ Xuefei Mi
}
\note{ No further notes

}
\seealso{
 No link
}
\examples{


Q=c(0.4308,0.9804,1.8603)
corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)



multistagegain(corr=corr,Q=Q)


#####
# code for using the parameter 'partial'
#####

multistagegain(corr=corr,Q=Q,Vg=0.4,partial=TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal }

