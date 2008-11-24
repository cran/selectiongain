\name{multistageselection}
\alias{multistageselection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ multistageselection}
\description{
 it is a function which calculates multi-stage selection gain.
 the gain can be calculated from a momentum function of multi-variate normaldistribution.
  }
\usage{
multistageselection(k, corr, alphaofx, sumdimofxandy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{ the lower bound of the integral, the first element of k is the lower bound of the response variable, Y, which is usually set as -Inf }
  \item{corr}{corrilation matrix }
  \item{alphaofx}{selected fraction }
  \item{sumdimofxandy}{ dimesion of k, which is dimesion of x plus 1(y has one dimension). }
}
\details{
 Tallis, 1961; Utz, 1984}
\value{
  the value returned, is the the gain of the selection.
  }
\references{ Tallis, 1965; Utz, 1984}
\author{ Xuefei, Mi}
\note{ further notes

}
\seealso{link}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(k,corr,alphaofx,sumdimofxandy)
{
 
dim=sumdimofxandy

A=array(0,c(dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
          if(i!=j)
         {
              A[i,j]= (k[j]-corr[i,j]*k[i])/ (1-corr[i,j]^2)^0.5
          }
    }
  }


part.corr=array(1,c(dim,dim,dim))

for (i in 1 : dim)
{  
   for (j in 1 : dim)
    {
         for (q in 1 : dim)
         {
         if(i!=j && q!=j && i!=q)

            { part.corr[i,j,q]= (corr[i,j]-corr[i,q]*corr[j,q])/ ((1-corr[i,q]^2)^0.5 * (1-corr[j,q]^2)^0.5)
             }
          }

      }
  }





j3q<-function (q,A,part.corr,dim)

{      
    
          
    lower=A[q,-q]
    corr= part.corr[-q,-q,q]


    output=pmvnorm(lower = lower, upper = rep(Inf,c(dim-1)), mean = rep(0, length(lower)), 
    corr = corr, sigma = NULL, algorithm = GenzBretz()) 
    output
 }

calculatx1<-function(A,part.corr,dim,corr,k,alpha3)
{ 
   output=0
   i=1
   for (i in 1 : dim)
 {
   output= output+ corr[1,i]*dnorm(k[i])*j3q(i,A,part.corr,dim)
 }
   
   output=output/alpha3

  output
  }
  

calculatx1(A=A,part.corr=part.corr,dim=dim,corr=corr,k=k,alpha3=alphaofx)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
