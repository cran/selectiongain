\name{multistageoptimal.pb}
\alias{multistageoptimal.pb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimizing two-stage selection in plant breeding with grid algorithm}
\description{This function is used to calculate the correlation matrix of n-stage selection for the testcross progenies of DH lines in all stages.
  
}
\usage{
multistageoptimal.pb(N.upper, N.lower, num.grid, Budget, 
CostC, CostTv, V, L1, Rep, N.fs, detail, fig, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 
\item{N.upper}{Vector with length n. It is the vector of upper limits of candidates X.}

\item{N.lower}{Vector with length n. It is the vector of lower limits of candidates X.}

\item{num.grid}{An integer value. It is the number of equally distanced points which divided the axis of \eqn{x_1} into \eqn{num.grid-1} intervals and there are \eqn{(num.grid-1)^n} grids in a n dimensional hyper cube.  More detail see \code{\link{multistageoptimal.grid}}.}

\item{Budget}{A double value. It contains the value of total budget.}

\item{CostC}{A double value. It contains the costs of producing and identifying a candidate.}

\item{CostTv}{Vector with length n. It contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}


  \item{V}{Vector of variance components. More detail see \code{\link{multistagecor}}.}  
 
  \item{L1}{Scalar of \eqn{L_1}.}
  
  \item{Rep}{Logical. If TRUE, the optimal selection index (Longin et al. 2007) will be used in the calculation. 
 }  
 

\item{N.fs}{An integer value. It is the number of final selected candidates.}

\item{detail}{Logical. if TRUE, the result of all the grids will be given or only the maximum. }

\item{alg}{An object used to switch between two algorithms. More detail see \code{\link{multistagegain}}.}

\item{fig}{Logical. If TRUE, a figure of contour plot will be saved in the default folder.}

}
\details{This is a special function made for plant breeding. A special and more complicated scenario relates to field tests of candidates in plant breeding. Here, the candidates are usually tested in replicated multi-location trials over several years, corresponding to the stages of selection. Thus, besides \eqn{\textbf{N}}, referring to the numbers of candidates to be tested in each stage, the breeder must also decide on the intensity of testing, as reflected by vector for the number of test locations \eqn{\textbf{L}=\{L_1,...,L_n\}} and replications \eqn{\textbf{R}=\{R_1,...,R_n\}}, where \eqn{L_i} and \eqn{R_i} refer to the number of test locations and replications per location, respectively, in stage \eqn{i}. 

If there is no upper limit on \eqn{L_i}, then \eqn{R_i}=1 is optimal for maximizing \eqn{\Delta G} (Longin et al. 2007). Normally a large number of candidates will be tested in few locations at the first stage, i.e., \eqn{L_1=1} or \eqn{2}. Under this scenario, the elements in \eqn{\bm{\Sigma}^{*}} are a rational function of \eqn{L}, \eqn{R} and the vector of variance components \eqn{\textbf{V}=\{Vg, Vgl, Vgy, Vgly, Ve \}}, where the latter refer to the variance among genotypes (\eqn{Vg}), genotype \eqn{\times} location interactions (\eqn{Vgl}), genotype \eqn{\times} year interactions (\eqn{Vgy}), genotype \eqn{\times} year \eqn{\times} location interactions (\eqn{Vgly}) and plot error (\eqn{Ve}). Here, \eqn{Vg} corresponds to \eqn{\sigma_y^2} and is set equal to 1. Likewise, the costs are not only a function of \eqn{\textbf{N}}, but also of \eqn{\textbf{L}} and \eqn{\textbf{R}}, because each test plot in field trials is associated with costs. Hence, the set of admissible allocations of resources \eqn{\Omega (B)}  can be described as


\deqn{
\Omega (B) :=\{ \omega = (\textbf{N,L,R})| C(\omega) \leq B \}.
}

In the simplest case,

\deqn{
 C(\omega) =  \sum_{i=1}^{n} N_i * L_i * R_i * CostT_i + N_1 * CostC \leq B,
}

where \eqn{CostC} refers to the costs of producing or identifying a candidate a candidate and \eqn{CostT_i} refers to the costs of testing a candidate in a test unit for stage \eqn{i}. }

\value{If detail = FALSE, the output of this function is a vector with the optimal number of candidates \eqn{(N_1,N_2)}, locations \eqn{(L_1,L_2)}, replicates \eqn{(R_1,R_2)} and the maximum \eqn{\Delta G(y)}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table.
}
%\references{

%}
\author{ Xuefei Mi}
\note{There is a limitation, if the absolute value of \eqn{q_i} is greater than 100, then the absolute value is settled to be 100. More detail about general application in optimization see \code{\link{multistageoptimal.grid}}.

}
\seealso{ \code{\link{multistagegain}}, \code{\link{multistageoptimal.grid}}}
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
  colnames(gain.table)= c("NumSelected","Budget","maxN1","maxN2",
  "Location1","steplength","Calcu.Gain")

gain.table[1,]=c(1,200, 101, 101,1,21 ,0)
gain.table[2,]=c(1,1000,401, 401,2,21,0)
gain.table[3,]=c(1,5000,2001,2001,2,41,0)
gain.table[4,]=c(4,200, 101, 101,1,21, 0)
gain.table[5,]=c(4,1000,601, 601,1,21, 0)
gain.table[6,]=c(4,5000,2001,2001,2,41,0) 

#######  
# IMPORTANT 
#######

# in order to reduce the time of checking in CRAN only the first breeding scheme will be checked

# if you want to run all 6 schemes you have to change the following code dim6=1 into dim6=6

dim6=1

#######  
# change the code above
#######


for (i in 1:dim6 )
{   
    maxn=gain.table[i,"maxN1"]
    length=gain.table[i,"steplength"]
    Budget=gain.table[i,"Budget"]
    Location1= gain.table[i,"Location1"]
    NumSelected=gain.table[i,"NumSelected"]    

    temp<-multistageoptimal.pb(N.upper=rep(maxn,2), N.lower=c(1,1), num.grid=length, 
    Budget=Budget,  CostC=0.5, CostTv=c(1,1), V="VC2", L1=Location1, Rep=c(1,1), 
    N.fs=NumSelected, alg=GenzBretz(),detail=TRUE)    
    result.grid[i,]=temp[[1]][1,]
}
  colnames(result.grid)<-c("NumSelected","Budget","Location1","Location2",
  "N1","N2","Rep1","Rep2","gain")
  rownames(result.grid)<-rep("grid",6)
  result.grid


# round-nlm

  dim=6
  gain.table= array(0,c(dim,7))
  result.nlm= array(0,c(dim,9))
  result.grid= array(0,c(dim,9))
  rownames(gain.table)= c(1:dim)
  colnames(gain.table)= c("NumSelected","Budget","maxN1","maxN2",
  "Location1","steplength","Calcu.Gain")

gain.table[1,]=c(1,200, 101, 101,1,21 ,0)
gain.table[2,]=c(1,1000,401, 401,2,21,0)
gain.table[3,]=c(1,5000,2001,2001,2,41,0)
gain.table[4,]=c(4,200, 101, 101,1,21, 0)
gain.table[5,]=c(4,1000,601, 601,1,21, 0)
gain.table[6,]=c(4,5000,2001,2001,2,41,0) 

# name changed in certain version

gainmatrix=gain.table

for (i in 1:dim6 )
{   
    maxn=gainmatrix[i,"maxN1"]
    length=gainmatrix[i,"steplength"]
    Budget=gainmatrix[i,"Budget"]
    Location1= gainmatrix[i,"Location1"]
    NumSelected=gainmatrix[i,"NumSelected"]    
    temp<-multistageoptimal.pb(N.upper=rep(maxn,2), N.lower=c(1,1), 
    num.grid=length, Budget=Budget,  CostC=0.5, CostTv=c(1,1), V="VC2", 
    L1=Location1, Rep=c(1,1), N.fs=NumSelected, alg=GenzBretz(),detail=FALSE)
    result.nlm[i,]=temp[[1]][2,]
}


  colnames(result.nlm)<-c("NumSelected","Budget","Location1","Location2",
  "N1","N2","Rep1","Rep2","gain")
  rownames(result.nlm)<-rep("round-nlm",6)
  result.nlm}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
