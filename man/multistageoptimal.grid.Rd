\name{multistageoptimal.grid}
\alias{multistageoptimal.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimizing with grid algorithm for fixed correlation matrix}
\description{This function calculates the maximum of \eqn{\Delta G(y)} with fixed correlation matrix by using grid algorithm.
  
}
\usage{
multistageoptimal.grid(N.upper, N.lower=rep(1,length(N.upper)), corr, 
num.grid=11, Budget, CostC=1, CostTv=rep(1,length(N.upper)),N.fs, 
detail=FALSE,alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 
\item{N.upper}{Vector with length n. It is the vector of uppper limits of candidates X.}

\item{N.lower}{Vector with length n. It is the vector of lower limits of candidates X.}


\item{corr}{(n+1-dimensional matrix). It is the correlation matrix \eqn{\bm{\Sigma}^{*}} of true value y and selection indices X. More detail see \code{\link{multistagegain}}.}

\item{num.grid}{ An integer value. It is the number of equally distanced points which divided the axis of \eqn{x_1} into \eqn{num.grid-1} intervals and there are \eqn{(num.grid-1)^n} grids in a n dimensional hyper cube. The default value of num.grid is 11, so the hyper cube is divided into \eqn{10^n} parts by the grid points. If the grid point is more than the maximum number of certain stage i, then the grid number of this stage will equal to the number of grid points, i.e., if \eqn{num.grid > N_i}, then the number of grid points for the i-th axis is \eqn{N_i}.}

\item{Budget}{A double value. It contains the value of total budget.}

\item{CostC}{A double value. It contains the costs of producing and identifying a candidate.}

\item{CostTv}{Vector with length n. It contains a vector with length n reflecting the cost of evaluating a candidate at stage i, i=1,...,n. The costs might vary in different stages.}

\item{N.fs}{An integer value. It is the number of final selected candidates.}

\item{detail}{ Logical. It is the control parameter to decide if the result of all the grids will be given (if TRUE) or only the maximum (if FALSE).}

\item{alg}{An object used to switch between two algorithms. More detail see \code{\link{multistagegain}}. }

}
\details{Suppose we start with \eqn{N_1} candidates in stage one. Based on the evaluation of the \eqn{N_i} candidates in stage \eqn{i}, the best \eqn{N_{i+1}} candidates, i.e., those with \eqn{x_i \geq q_i}, are promoted to the next stage, where they are evaluated with even higher intensity to obtain more precise estimates of the true value \eqn{y}. The goal of the whole selection scheme is to select the best \eqn{N_{n+1}} candidates after \eqn{n} stages of selection.

In practice, the selection program has only a limited budget \eqn{B} to cover all costs such as (i) identifying or producing the initial \eqn{N_1} candidates and (ii) evaluating the \eqn{N_i} candidates in stage \eqn{i}. For a given testing scheme with \eqn{\textbf{N}=(N_1,\ldots,N_n)} candidates in the \eqn{i}-th stage of selection (\eqn{i=1,\ldots,n}), the costs may be given by the cost function \eqn{C(\omega)}. Thus, the set of admissible allocations \eqn{\Omega (B)} of the candidates to the various stages of selection is given by 


\deqn{
\Omega (B):= \{ \omega =\textbf{N}|C(\omega)\leq B\}
} 
 
Hence, our goal is to find \eqn{\tilde{\omega} \in \Omega (B)} with 

\deqn{
\Delta G(y|\textbf{S}_{\tilde{\omega}}, \Sigma^{*}) = \underbrace{max}_{\omega \in \Omega (B)} \Delta G(y|\textbf{S}_{\omega},\Sigma^{*}),
} 
 
where \eqn{\textbf{S}_{\omega}} refers to the truncation point \eqn{\textbf{Q}} corresponding to \eqn{\textbf{A}=\{ \alpha_1,\ldots,\alpha_n \}}, with \eqn{\alpha_i=N_{i+1}/N_i} for \eqn{i=1,\ldots,n}. The matrix \eqn{ \bm{\Sigma}^{*}} is determined by the correlations among test scores \eqn{x_i} obtained in the \eqn{n} stages of selection as well as their correlations to the target value \eqn{y}. Hence, for given but possibly different testing procedures in each stage,  \eqn{\bm{\Sigma}^{*}} is fixed, independent of the choice of \eqn{\textbf{N}}.


In many applications in breeding and other fields, the choice of \eqn{\textbf{N}} does not affect the correlation matrix \eqn{\bm{\Sigma}^{*}} for the candidates. Examples include different types of average in the various stages of selection such as tests for evaluating the disease symptoms (e.g., test of fusarium resistance by visual recording of disease symptoms, estimation of mycotoxin concentration by NIRS, ELISA or GC-MS) or genomic usages with different prediction accuracy and costs (marker arrays with different coverage of the genome, transcript and for metabolic profiles). All these situations can be coped within this frame work outlined above.


The simplest way to find the maximum is to do a full scan of the entire set \eqn{\Omega (B)}, which calculates \eqn{\Delta G(y|\textbf{S}_{\omega}, \bm{\Sigma}^{*}) } for all possible allocations of \eqn{\omega (B)} to determine \eqn{\tilde{\omega}} yielding the largest \eqn{\Delta G}. However, this is very time consuming. An alternative solution is to use grid search, which divides the whole set \eqn{\Omega (B)} into several grids (Kim 1997).

}
\value{
If \eqn{\texttt{detail}} = FALSE, the output of this functions is a vector with the optimal number of candidates in each stage (\eqn{\textbf{N}}) and the maximum \eqn{\Delta G(y)}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table.}
\references{
Kim, J. (1997). Iterated Grid Search Algorithm on Unimodal Criteria. Ph.D. thesis, Virginia Polytechnic Institute and State University.
}
\author{ Xuefei Mi}
\note{No further comment

}
\seealso{ \code{\link{multistagegain}} }
\examples{

 
corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)

multistageoptimal.grid(N.upper=rep(100,3), corr=corr, Budget=200, CostC=0.5, N.fs=5)


correlation=matrix(c(
   1,0.2,0.3,0.5,
   0.2,1,0,0,
   0.3,0,1,0,
   0.5,0,0,1 ),
4,4)

CostV=c(1,10,20)

multistageoptimal.grid(N.upper=c(1000,200,100)+1, corr=correlation,
Budget=2200, CostC=0, CostTv=CostV, N.fs=5)


#######  
# IMPORTANT 
#######

# in order to reduce the time of checking in CRAN only the first breeding scheme will be checked

# if you want to run all 5 schemes you have to change the following code, 
# just simply remove the number sign.

#######  
# change the code below
#######


#multistageoptimal.nlm(N.upper=c(1000,200,100), corr=correlation, 
#Budget=2200, CostC=0, CostTv=CostV, N.fs=5)

#multistageoptimal.nlm(N.upper=c(1000,200,100), corr=correlation, 
#Budget=2200, CostC=0, CostTv=CostV, N.fs=5,ini.value=c(701,81,31))

#multistageoptimal.grid(N.upper=c(1000,200,100)+1, corr=correlation, 
#Budget=2200, CostC=0, CostTv=CostV, N.fs=5,num.grid=21)

#multistageoptimal.nlm(N.upper=c(1000,200,100), corr=correlation, 
#Budget=2200, CostC=0, CostTv=CostV, N.fs=5,ini.value=c(800,100,20))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
