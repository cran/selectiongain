\name{multistagevariance}
\alias{multistagevariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected gain for k-stages selection}
\description{
This function uses the algorithm described by Tallis (1961) to calculate the variance (second moment) after multi-stage selection. 

}
\usage{
multistagevariance(Q, corr, alg, lim.y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
\item{Q}{ Vector with length n. It refers to the coordinates of the truncation points Q, which is the output of the function \code{\link{multistagetp}} that we are going to introduced.}
  
\item{corr}{(n+1)-dimensional matrix. It is the correlation matrix of true value y and selection indices X.  More detail see \code{\link{multistagegain}}.}

\item{alg}{An object used to switch between two algorithms. More detail see \code{\link{multistagegain}}.}
 
\item{lim.y}{The lower limit of y as double, and set to -200 as default, which is on the safe side.} 
   
}


\details{

The theory for multi-stage selection is based on the first moment (mathematical expectation) of the true value \eqn{y} in the selected area Cochran (1951). The \eqn{\Delta G(y)} is defined as the difference of the mathematical expectation of \eqn{y} after selection and the corresponding value before selection.

  For simplicity, the variance of \eqn{y} is set to \eqn{1}, i.e., \eqn{\sigma_{y}^2=1}. If \eqn{\sigma_{y}^2 \neq1}, the selection gain and the variance among the selected candidates have to be multiplied with \eqn{\sigma_{y}} or \eqn{\sigma_{y}^2}, respectively. In order to calculate the expectation of \eqn{y} after \eqn{n} stages of selection, we have to determine the one-sided integral of \eqn{y} over the right-sided area \eqn{\textbf{S}_{Q}=\{x_1>q_1,\ldots,x_n>q_n\}} defined by the truncation point \eqn{\textbf{Q}=\{q_1,q_2,\ldots,q_n\}}. The value of \eqn{\Delta G_n(y)} is noted as \eqn{E(Y;\textbf{S}_{Q})}. We have to distinguish it with \eqn{E(Y)}, which is the integral of \eqn{y} over the whole area \eqn{\textbf{S}=\{x_1>-\infty,\ldots,x_n>-\infty\}} with \eqn{\textbf{Q}=\{-\infty,\ldots,-\infty\}} as the truncation point. Thus, the expectation of \eqn{y} after \eqn{n} stages of selection is,

\deqn{
 \Delta G_n(y, \textbf{S}_{Q}) = E(Y;\textbf{S}_{Q})= \alpha^{-1} \int_{-\infty} ^\infty   \int_{q_{1}}^\infty\ldots\int_{q_{n}}^\infty y \,  \phi_{n+1}(\textbf{x}^{*};\textbf{U}^{*}; \bm{\Sigma}^{*}) \, d \textbf{x}^*, 
}

where
\deqn{ 
\alpha = \Phi_n (\textbf{Q},\bm{\Sigma})= \int_{q_{1}}^\infty\ldots\int_{q_{n}}^\infty \phi_{n}(\textbf{x};\textbf{U}; \bm{\Sigma}) \, d \textbf{x}, \label{equ-alpha0}
}


and \eqn{\Phi_n} is the distribution function of the MVN, \eqn{\phi_{n}} is the density function of MVN, and \eqn{\bm{\Sigma}} is the correlation matrix of \eqn{\textbf{X}}. \eqn{\bm{\Sigma}^{*}} is the correlation matrix of \eqn{\textbf{X}^*}. It comprises \eqn{\bm{\Sigma}}, but has one dimension more pertaining to the correlations between \eqn{Y=X_0} and the selection indices \eqn{\textbf{X}}.  The mean vector \eqn{\textbf{U}^{*}=\{u_0,u_1,\ldots,u_n\}} of \eqn{\phi_{n+1}} is omitted, assuming \eqn{\textbf{U}^{*}=\{0,0,\ldots,0\}} without loss of generality, and consequently, we write shortly \eqn{\phi_{n+1} (\textbf{x}^{*};\bm{\Sigma}^{*})} and \eqn{\phi_{n} (\textbf{x};\bm{\Sigma})}. The selection gain is the first moment, while the selected fraction \eqn{\alpha} over all \eqn{n} stages of selection corresponds to the zero-th moment of the one-sided truncated MVN distribution of \eqn{\textbf{X}}. 

The variance among candidates of y in the selected area \eqn{\textbf{S}_{Q}} is defined as the second central moment, \eqn{\psi_n(y)=E(Y^2;\textbf{S}_{Q}) - \left[E(Y;\textbf{S}_{Q})\right]^2}, where

\deqn{
E(Y^2;\textbf{S}_{Q})  = \alpha^{-1} \int_{-\infty} ^\infty \int_{q_{1}}^\infty...\int_{q_{n}}^\infty  y^2\, \phi_{n+1}(\textbf{x}^{*}; \bm{\Sigma}^{*}) \, d \textbf{x}^*}

}

\value{
The output is the value of \eqn{\psi_n(y;\textbf{S}_{Q})=E(y^2;\textbf{S}_{Q})-[E(y;\textbf{S}_{Q})]^2}.
  }
%\references{}

\author{Xuefei Mi
}
\note{No further notes

}
\seealso{
 \code{\link{multistagegain}}
}
\examples{
Q =c(0.4308,0.9804,1.8603)

corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)


multistagevariance(Q=Q,corr=corr,alg=Miwa)

var.time.miwa=system.time (var.miwa<-multistagevariance(Q=Q,corr=corr,alg=Miwa))

var.time.bretz=system.time (var.bretz<-multistagevariance(Q=Q,corr=corr))


var.time.miwa
var.miwa[1]
var.time.bretz
var.bretz[1]

# further examples 1 


Q= c(0.9674216, 1.6185430)
corr=matrix( c(1,      0.7071068, 0.9354143,
               0.7071068, 1,      0.7559289,
               0.9354143, 0.7559289, 1    
             ), 
              nrow=3  
)


multistagevariance(Q=Q,corr=corr,alg=Miwa)

var.time.miwa=system.time (var.miwa<-multistagevariance(Q=Q,corr=corr,alg=Miwa))

var.time.bretz=system.time (var.bretz<-multistagevariance(Q=Q,corr=corr))


var.time.miwa
var.miwa[1]
var.time.bretz
var.bretz[1]


# further examples 

 alpha1<- 1/(24)^0.5
 alpha2<- 1/(24)^0.5
 Q=multistagetp(alpha=c(alpha1,alpha2),corx=corr[2:3,2:3])


corr=matrix( c(1,      0.7071068, 0.9354143,
               0.7071068, 1,      0.7559289,
               0.9354143, 0.7559289, 1    
             ), 
              nrow=3  
)

multistagevariance(Q=Q,corr=corr,alg=Miwa)






}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ Truncated multivariate normal }% __ONLY ONE__ keyword per line
