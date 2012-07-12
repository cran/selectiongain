\name{multistageoptimal.pb}
\alias{multistageoptimal.pb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for optimizing two-stage selection in plant breeding with a
8 grid algorithm}
\description{This function is used to calculate the correlation matrix of $n$-stage selection for the testcross progenies of DH lines in all stages.
  
}
\usage{
multistageoptimal.pb(N.upper, N.lower, num.grid, Budget, 
CostC, CostTv, V, L1, Rep, N.fs, detail, fig, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 
\item{N.upper}{ is the vector of up limits of number of candidates X.}

\item{N.lower}{is the vector of low limits of number of candidates X.}

\item{num.grid}{ is the number of equally distanced points which divided the axis of \eqn{x_1} into \eqn{num.grid-1} intervals and there are \eqn{(num.grid-1)^n} grids in a n dimensional hyper cube. The default value of num.grid is 11, so the hyper cube is divided into \eqn{10^n} parts by the grid points. If the grid point is more than the maximum number of certain stage i, then the grid number of this stage will equal to the number of grid points, i.e., if \eqn{num.grid > N_i}, then the number of grid points for the i-th axis is \eqn{N_i}.}


\item{Budget}{contains the value of total budget.}

\item{CostC}{ contains the costs of producing or identifying a candidate.}

\item{CostTv}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}


  \item{V}{is the vector of variance components.}  
 
  \item{L1}{is the scalar of \eqn{L_1}.}
  
  \item{Rep}{is the control parameter. If it equals TRUE, the optimal selection index (longin2007) will be used in the calculation. 
 }  
 

\item{N.fs}{is the number of final selected candidates.}

\item{detail}{is the control parameter to decide if the result of all the grids will be given or only the maximum. }

\item{alg}{is used to switch between two algorithms. If alg = GenzBretz(), which is by default, the quasi-Monte Carlo algorithm from Genz(1999) will be used. If alg = Miwa(), the program will use the Miwa algorithm (Mi2009), which an analytical solution of the MVN integral (Miwa2003). Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits), however the computational speed is also slower. Therefore, we recommend the user to use the default value of this parameter.}

\item{fig}{is the control parameter to decide if a figure of contour plot will be saved in the default folder.}

}
\details{There is a limitation, if the absolute value of \eqn{q_i} is greater than 100, then the absolute value is settled to be 100.}
\value{If detail = FALSE, the output of this function is a vector with the optimal number of candidates \eqn{(N_1,N_2)}, locations \eqn{(L_1,L_2)}, replicates \eqn{(R_1,R_2)} and the maximum \eqn{\Delta G(y)}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table.
}
\references{
W.G. Cochran. Improvent by means of selection. In: Proceedings Second Berkeley Symposium on Math Stat Prof, pp449-470.

G.M. Tallis. Moment generating function of truncated multi-normal distribution. Journal of the Royal Statistical Society, Series B, 23(1):223-229, 1961.

H.F. Utz. Mehrstufenselecktion in der Pflanzenzuechtung. Doctor thesis, University Hohenheim, 1969.

X. Mi, H.F. Utz. and A.E. Melchinger.  R package selectiongain: A tool for efficient calculation and optimization of the expected gain from multi-stage selection. J Stat Softw. (submitted)
}
\author{ Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{

# examples for the JSS paper

multistageoptimal.pb(N.upper=rep(401,2), N.lower=c(1,1), num.grid=21, Budget=1000, 
 CostC=0.5, CostTv=c(1,1), V="VC2", L1=2, Rep=c(1,1), N.fs=1, alg=GenzBretz())

# glm

  dim=6
  gain.table= array(0,c(dim,7))
  result.nlm= array(0,c(dim,9))
  result.grid= array(0,c(dim,9))
  rownames(gain.table)= c(1:dim)
  colnames(gain.table)= c("NumSelected","Budget","maxN1","maxN2","Location1","steplength","Calcu.Gain")

gain.table[1,]=c(1,200, 101, 101,1,21 ,0)
gain.table[2,]=c(1,1000,401, 401,2,21,0)
gain.table[3,]=c(1,5000,2001,2001,2,41,0)
gain.table[4,]=c(4,200, 101, 101,1,21, 0)
gain.table[5,]=c(4,1000,601, 601,1,21, 0)
gain.table[6,]=c(4,5000,2001,2001,2,41,0) 


for (i in 1:dim )
{   
    maxn=gain.table[i,"maxN1"]
    length=gain.table[i,"steplength"]
    Budget=gain.table[i,"Budget"]
    Location1= gain.table[i,"Location1"]
    NumSelected=gain.table[i,"NumSelected"]    

    temp<-multistageoptimal.pb(N.upper=rep(maxn,2), N.lower=c(1,1), num.grid=length, Budget=Budget,  CostC=0.5, CostTv=c(1,1), V="VC2", L1=Location1, Rep=c(1,1), N.fs=NumSelected, alg=GenzBretz(),detail=TRUE)    
    result.grid[i,]=temp[[1]][1,]
}
  colnames(result.grid)<-c("NumSelected","Budget","Location1","Location2","N1","N2","Rep1","Rep2","gain")
  rownames(result.grid)<-rep("grid",6)
  result.grid


# round-nlm

  dim=6
  gain.table= array(0,c(dim,7))
  result.nlm= array(0,c(dim,9))
  result.grid= array(0,c(dim,9))
  rownames(gain.table)= c(1:dim)
  colnames(gain.table)= c("NumSelected","Budget","maxN1","maxN2","Location1","steplength","Calcu.Gain")

gain.table[1,]=c(1,200, 101, 101,1,21 ,0)
gain.table[2,]=c(1,1000,401, 401,2,21,0)
gain.table[3,]=c(1,5000,2001,2001,2,41,0)
gain.table[4,]=c(4,200, 101, 101,1,21, 0)
gain.table[5,]=c(4,1000,601, 601,1,21, 0)
gain.table[6,]=c(4,5000,2001,2001,2,41,0) 

# name changed in certain version

gainmatrix=gain.table

for (i in 1:dim )
{   
    maxn=gainmatrix[i,"maxN1"]
    length=gainmatrix[i,"steplength"]
    Budget=gainmatrix[i,"Budget"]
    Location1= gainmatrix[i,"Location1"]
    NumSelected=gainmatrix[i,"NumSelected"]    
    temp<-multistageoptimal.pb(N.upper=rep(maxn,2), N.lower=c(1,1), num.grid=length, Budget=Budget,  CostC=0.5, CostTv=c(1,1), V="VC2", L1=Location1, Rep=c(1,1), N.fs=NumSelected, alg=GenzBretz(),detail=FALSE)
    result.nlm[i,]=temp[[1]][2,]
}


  colnames(result.nlm)<-c("NumSelected","Budget","Location1","Location2","N1","N2","Rep1","Rep2","gain")
  rownames(result.nlm)<-rep("round-nlm",6)
  result.nlm}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
