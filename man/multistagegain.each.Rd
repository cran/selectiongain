\name{multistagegain.each}
\alias{multistagegain.each}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the selection gain in each stage}
\description{
In some situations, the user wants to know the increase of \eqn{\Delta G} for each stage so that it is possible to determine the stage which contributes most to \eqn{\Delta G}. This function calculates \eqn{\Delta G} stepwise for each stage. 
}
\usage{
multistagegain.each(corr, Q, Vg, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{


\item{corr}{is the correlation matrix of y and X, which is introduced in function multistagecorr. The correlation matrix must be symmetric and positive-definite. Before starting the calculations, the user is recommended to check the correlation matrix.}

\item{Q}{are the coordinates of the truncation points, which are the output of the next function multistagetp that we are going to introduce.}
  
\item{Vg}{is true breeding value. The default value is 1.}


\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz (2009) will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et. al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, the computational speed is slower. We recommend the user to use Miwa algorithm of this parameter.}
 
 
}
\details{This function calculates the well-known selection gain \eqn{\Delta G}, which is described by Cochran (1951). For one-stage selection the gain is defined as \eqn{\Delta G = i V_g \rho_{1}}, where \eqn{i} is the selection intensity, \eqn{\rho_{1}} is the correlation between the true breeding value (its variance is \eqn{V_g}) and the selection index (Utz 1969).

}

\value{
  The output is given as \eqn{(\Delta G_1(y), \Delta G_2(y)-\Delta G_1(y), \Delta G_3(y)-\Delta G_2(y), ...)$ where $\Delta G_i(y)} refers to the total selection gain after the first i stages of selection. 
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
%\note{ No further notes

%}
\seealso{
selectiongain()
}
\examples{


corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)

multistagegain.each(Q=c(0.4308,0.9804,1.8603),corr=corr)

# further examples 3 for the JSS paper

 alpha1<- 1/24
 alpha2<- 1
 Q=multistagetp(alpha=c(alpha1,alpha2),corx=corr[2:3,2:3])

k=c(-200,Q)

corr=matrix( c(1,      0.7071068, 0.9354143,
               0.7071068, 1,      0.7559289,
               0.9354143, 0.7559289, 1    
             ), 
              nrow=3  
)

alphaofx=pmvnorm(lower=k,corr=corr)

multistagegain(Q=Q,corr=corr)

multistagegain(Q=Q,corr=corr,partial=TRUE)

multistagegain.each(Q=Q,corr=corr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal }

