\name{multistagegain}
\alias{multistagegain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the multi-stage selection gain}
\description{
 
This is the main function of the package and calculates selection gain for a fixed correlation matrix and coordinates of the truncation points. 

}
\usage{
multistagegain(Q, corr, alg, lim.y, stages=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 
\item{Q}{ Vector with length n. It refers to the coordinates of the truncation points Q, which is the output of the function \code{\link{multistagetp}}.}
  
\item{corr}{(n+1)-dimensional matrix. It is the correlation matrix \eqn{\bm{\Sigma}^{*}} of true value y and selection indices X. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix, which is usually obtained by analysis of one or several experiments designed for estimating covariance components as well as phenotypic and genotypic correlations.}

\item{alg}{ An object used to switch between two algorithms. If alg = GenzBretz(), which is by default, the quasi-Monte Carlo algorithm from Genz (1999) will be used. If alg = Miwa(), the program will use the Miwa algorithm, which an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits), however the computational speed is also slower. Therefore, we recommend the user to use the default value of this parameter.}
 
\item{lim.y}{ The lower limit of y as double, and set to -200 as default, which is on the safe side.} 
  
\item{stages}{Logical. If TRUE, the Tallis' partial gain of deleting one stage from n stages, which are integrals with partial correlation matrix, will be shown (Tallis 1961). The sum of these partial gains will be equal to the selection gain. By default it equals FALSE.}



}
\details{

This function calculates the selection gain \eqn{\Delta G(y)}, which is described by Cochran (1951). For one-stage selection the gain is defined as \eqn{\Delta G(y) = i \rho_{y} \rho_{1}}, where \eqn{i} is the selection intensity, \eqn{\rho_{1}} is the correlation between the true breeding value and the selection index \eqn{y} (Utz 1969). The theory for multi-stage selection is based on the first moment (mathematical expectation) of the true value \eqn{y} in the selected area Cochran (1951). The \eqn{\Delta G(y)} is defined as the difference of the mathematical expectation of \eqn{y} after selection and the corresponding value before selection.
 
  For simplicity, the variance of \eqn{y} is set to \eqn{1}, i.e., \eqn{\sigma_{y}^2=1}. If \eqn{\sigma_{y}^2 \neq1}, the selection gain and the variance among the selected candidates have to be multiplied with \eqn{\sigma_{y}} or \eqn{\sigma_{y}^2}, respectively. In order to calculate the expectation of \eqn{y} after \eqn{n} stages of selection, we have to determine the one-sided integral of \eqn{y} over the right-sided area \eqn{\textbf{S}_{Q}=\{x_1>q_1,\ldots,x_n>q_n\}} defined by the truncation point \eqn{\textbf{Q}=\{q_1,q_2,\ldots,q_n\}}. The value of \eqn{\Delta G_n(y)} is noted as \eqn{E(Y;\textbf{S}_{Q})}. We have to distinguish it with \eqn{E(Y)=0}, which is the integral of \eqn{y} over the whole area \eqn{\textbf{S}=\{x_1>-\infty,\ldots,x_n>-\infty\}} with \eqn{\textbf{Q}=\{-\infty,\ldots,-\infty\}} as the truncation point. Thus, the expectation of \eqn{y} after \eqn{n} stages of selection is,

\deqn{
 \Delta G_n(y, \textbf{S}_{Q}) = E(Y;\textbf{S}_{Q})= \alpha^{-1} \int_{-\infty} ^\infty   \int_{q_{1}}^\infty\ldots\int_{q_{n}}^\infty y \,  \phi_{n+1}(\textbf{x}^{*};\textbf{U}^{*}; \bm{\Sigma}^{*}) \, d \textbf{x}^*, 
}

where
\deqn{ 
\alpha = \Phi_n (\textbf{Q},\bm{\Sigma})= \int_{q_{1}}^\infty\ldots\int_{q_{n}}^\infty \phi_{n}(\textbf{x};\textbf{U}; \bm{\Sigma}) \, d \textbf{x}, \label{equ-alpha0}
}


and \eqn{\Phi_n} is the distribution function of the MVN, \eqn{\phi_{n}} is the density function of MVN, and \eqn{\bm{\Sigma}} is the correlation matrix of \eqn{\textbf{X}}. \eqn{\bm{\Sigma}^{*}} is the correlation matrix of \eqn{\textbf{X}^*=\{Y, X\}}. It comprises \eqn{\bm{\Sigma}}, but has one dimension more pertaining to the correlations between \eqn{Y=X_0} and the selection indices \eqn{\textbf{X}}.  The mean vector \eqn{\textbf{U}^{*}=\{u_0,u_1,\ldots,u_n\}} of \eqn{\phi_{n+1}} is omitted, assuming \eqn{\textbf{U}^{*}=\{0,0,\ldots,0\}} without loss of generality, and consequently, we write shortly \eqn{\phi_{n+1} (\textbf{x}^{*};\bm{\Sigma}^{*})} and \eqn{\phi_{n} (\textbf{x};\bm{\Sigma})}. The selection gain is the first moment, while the selected fraction \eqn{\alpha} over all \eqn{n} stages of selection corresponds to the zero-th moment of the one-sided truncated MVN distribution of \eqn{\textbf{X}}. 



}

\value{
  The value returned, is the gain of selection.
  }
\references{

Cochran, W.G. (1951). Improvement by Means of Selection. \emph{Proc. Second Berkeley Symp. on Math. Statist. and Prob.}, 449-470. University of California Press, Berkeley.

Genz, A., Bretz F., Miwa, T., Mi, X., Leisch, F., Scheipl, F. and Hothorn, T. (2012) mvtnorm: Multivariate normal and t distributions, R package version 0.9-9994.

Tallis, G.M. (1961) Moment generating function of truncated multi-normal distribution. \emph{Journal of the Royal Statistical Society, Series B},  \bold{23(1)}, 223-229.

Mi, X., Miwa, T. and Hothorn T. (2009) Implement of Miwa's analytical algorithm of multi-normal distribution, \emph{R Journal}, \bold{1}, 37-39.

Mi, X., Utz, H.F. and Melchinger, A.E. (submitted). R Package selectiongain: Calculating and Optimizing Expected Gain from Multi-stage Selection. \emph{Journal of Statistical Software} 


}


\author{ Xuefei Mi
}
\note{ 
Equation given by Tallis (1961) are used for calculating the above integral. The partial gains are gains if one stage is deleted, are not equal to the gain from each stages, which is calculated by \code{\link{multistagegain.each}}.

}
\seealso{\code{\link{pmvnorm}}
}
\examples{


Q=c(0.4308,0.9804,1.8603)

corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)

multistagegain(Q=Q,corr=corr)

# code for partial gain

multistagegain(Q=Q,corr=corr,stages=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal }

