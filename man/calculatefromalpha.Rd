\name{calculatefromalpha}
\alias{calculatefromalpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calculate the quantile for given alpha}
\description{
  this function can be considered as similar as qmvnorm. but for qmvnorm, a uniform alpha is given and the uniform quantiles are calculated for this alpha.
in our function a vector of alpha is given, and the quantiles are calculate one by one to satisfy thesis given alpha. 

 e.g.  alpha2= pmvnorm(quantile1,quantile2)  
    now given alpha2, quantile1. we want to know quantile2.
    this is kind of Recursion Procedure.

}
\usage{
calculatefromalpha(alpha, dim, corr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{  Probability vector for random variable X. $alpha(i)= pmvnorm(quantile(1),..., quantile(i)) $
        in plants breeding, it is also called selection fraction.  here we have $alpha(i)=N_(i+1)/N(i)$,  which is the elimation rate from generation $i$ to genereation $i+1$.
 }
  \item{dim}{ dimesion of the random variable X.}
  \item{corr}{ correlation matrix }
}
\details{
  more details
}
\value{
  the  out come is a vector for the quantiles
}
\references{ Tallis 1961, Utz1969}
\author{ Xuefei Mi}
\note{ further notes 

}
\seealso{ selectiongain() }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(alpha,dim,corr)
{

    qnormbar<-function(kwanted,koutput,corr,howmanyk,alpha)
   {   
        corr2dim=corr[1:howmanyk,1:howmanyk] 
        
        tempofk1=koutput[1:c(howmanyk-1)]
          
        alphatotal=1
        for (i in 1:howmanyk)
        {
           alphatotal=alphatotal*alpha[i]
        }

        pmvnorm(lower = c(tempofk1,kwanted), upper = rep(Inf,howmanyk), mean = rep(0,howmanyk),
                   corr = corr2dim) - alphatotal   
   }
     
    koutput=0
    for (i in  1: dim )
       {
    
         if (i == 1)     
         {
           koutput[1]= qnorm(alpha[1],lower.tail = FALSE)[[1]][1]
  
         }else
         {

           koutput[i]= uniroot(qnormbar,interval=c(-100,100),koutput=koutput,corr,i,alpha)[[1]][1]  
         }
     }

koutput
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
