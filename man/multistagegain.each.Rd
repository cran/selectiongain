\name{multistagegain.each}
\alias{multistagegain.each}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating the selection gain in each stage}
\description{

 This function calculates \eqn{\Delta G(y)} stepwise for each stage. 
}
\usage{
multistagegain.each(Q, corr, alg, lim.y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{


\item{Q}{ Vector with length n. It refers to the coordinates of the truncation points Q, which is the output of the function \code{\link{multistagetp}}.}
  
\item{corr}{(n+1-dimensional matrix). It is the correlation matrix \eqn{\bm{\Sigma}^{*}} of true value y and selection indices X. More detail see \code{\link{multistagegain.each}}.}


\item{alg}{An object used to switch between two algorithms. More detail see \code{\link{multistagegain}}.}
 
\item{lim.y}{The lower limit of y as double, and set to -200 as default, which is on the safe side.} 

}
\details{In some situations, the user wants to know the increase of \eqn{\Delta G(y)} in each stage so that it is possible to determine the stage which contributes most to \eqn{\Delta G(y)}.
}

\value{
  The output is given as \eqn{(\Delta G_1(y), \Delta G_2(y)-\Delta G_1(y), \Delta G_3(y)-\Delta G_2(y), ...)} where \eqn{\Delta G_i(y)} refers to the total selection gain after the first i stages of selection. 
  }
%\references{

%}


\author{ Xuefei Mi
}
\note{ No further notes

}
\seealso{
\code{\link{multistagegain}}
}
\examples{
k=c(-200,0.4308,0.9804,1.8603)

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

multistagegain(Q=Q,corr=corr,)

multistagegain(Q=Q,corr=corr,stages=TRUE)

multistagegain.each(Q=Q,corr=corr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated multivariate normal }

