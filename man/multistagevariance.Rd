\name{multistagevariance}
\alias{multistagevariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected variance after selection for k-stages selection}
\description{
This function uses the algorithm described by Tallis (1961) to calculate the variance after multi-stage selection. 
The variance among candidates of y in the selected area \eqn{\textbf{S}_{Q}} is defined as the second central moment, \eqn{\psi_n(y)=E(Y^2|\textbf{S}_{Q}) - \left[E(Y|\textbf{S}_{Q})\right]^2}, where
\eqn{
E(Y^2|\textbf{S}_{Q})  = \alpha^{-1} \int_{-\infty} ^\infty \int_{q_{1}}^\infty...\int_{q_{n}}^\infty  y^2\, \phi_{n+1}(\textbf{x}^{*}; \bm{\Sigma}^{*}) \, d \textbf{x}^*}


}
\usage{
multistagevariance(Q, corr, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Q}{ (length n) is the truncation points Q, which is the output of the next function multistagetp that we are going to introduce.}
  
  \item{corr}{ is the correlation matrix of y and X, which is introduced in function multistagecorr. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix. }

\item{alg}{is used to switch between two algorithms. If alg = GenzBretz(), which is by default, the quasi-Monte Carlo algorithm from Genz (1999) will be used. If alg = Miwa(), the program will use the Miwa algorithm (Mi et. al., 2009), which an analytical solution of the MVN integral (Miwa, 2003). Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits), however the computational speed is also slower. Therefore, we recommend the user to use the default value of this parameter.}

}
%\details{
%More details are in the JSS paper section 3.4.
%}

\value{
The output is the value of \eqn{\psi_n(y|\textbf{S}_{Q})}.
}
  
  
  
\references{

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution, R Journal, 1:37-39, 2009.

X. Mi, H.F. Utz. and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)
}


\author{Xuefei Mi
}
\note{ No further notes

}
\seealso{
 No link
}
\examples{

# first example

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

# second examples 


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


# third examples 

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
