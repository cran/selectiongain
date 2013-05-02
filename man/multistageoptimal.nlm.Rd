\name{multistageoptimal.nlm}
\alias{multistageoptimal.nlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing n-stage selection with the NLM algorithm for fixed correlation matrix}
\description{This function uses the non-linear minimization function nlm in R-package stat for n-stage selection
  }
\usage{
multistageoptimal.nlm(N.upper, N.lower, corr, ini.value, 
Budget, CostC, CostTv, N.fs, iterlim, alg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{N.upper}{Vector with length n. It is the vector of upper limits of candidates X.}

\item{N.lower}{Vector with length n. It is the vector of lower limits of candidates X.}


\item{corr}{(n+1)-dimensional matrix. It is the correlation matrix \eqn{\bm{\Sigma}^{*}} of true value y and selection indices X. More detail see \code{\link{multistagegain.each}}.}

\item{ini.value}{Vector with length n. It stores the number of candidates in each stage for the algorithm to begin with. As default, it will use \eqn{N=\{N_1,N_2,...,N_n\}=\{a+1,...,a+n\}}, where a is defined as \code{Budget/(CostC+sum(CostTv)+1).}}


\item{Budget}{A double value. It contains the value of total budget.}

\item{CostC}{A double value. It contains the costs of producing and identifying a candidate.}

\item{CostTv}{A double value. It contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages.}

\item{N.fs}{Vector with length n. It is the number of final selected candidates.}


\item{iterlim}{An integer value. It is the maximum number of iterations to be executed before the Newton algorithm is terminated. By default it is equal to 20. If the \eqn{\texttt{Budget}} increases 10 times for making the selection, the value of \eqn{\texttt{iterlim}} has to be increased \eqn{lg(10)} times.}

\item{alg}{An object used to switch between two algorithms. More detail see \code{\link{multistagegain}}.}
}
\details{

Suppose we start with \eqn{N_1} candidates in stage one. Based on the evaluation of the \eqn{N_i} candidates in stage \eqn{i}, the best \eqn{N_{i+1}} candidates, i.e., those with \eqn{x_i \geq q_i}, are promoted to the next stage, where they are evaluated with even higher intensity to obtain more precise estimates of the true value \eqn{y}. The goal of the whole selection scheme is to select the best \eqn{N_{n+1}} candidates after \eqn{n} stages of selection.



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


The simplest way to find the maximum is to do a full scan of the entire set \eqn{\Omega (B)}, which calculates \eqn{\Delta G(y|\textbf{S}_{\omega}, \bm{\Sigma}^{*}) } for all possible allocations of \eqn{\omega (B)} to determine \eqn{\tilde{\omega}} yielding the largest \eqn{\Delta G}. However, this is very time consuming. An alternative solution is to use grid search, which divides the whole set \eqn{\Omega (B)} into several grids Kim (1997).
 
  Another way for finding the maximum is using an optimization algorithm for Non-Linear Minimization (NLM) provided by  function \code{\link{nlm}} in package \code{\link{stats}}. It uses a Newton-type algorithm for searching the maximum of a multi-modal function. This algorithm depends heavily on the starting point, the maximum number of iterations and the numerical derivatives of  \eqn{\Delta G} and results in an accuracy less than four digits. For maximizing the selection gain, NLM algorithm will converge to the global maximum . A proper choice of the initial value is recommended.

}
\value{
The output of this function is a vector similar as in \code{\link{multistageoptimal.grid}} if detail = FALSE. However, the optimal number of candidates in each stage determined by the NLM algorithm is clearly not an integer, because the function uses a numerical algorithm, which depends on derivatives. 
}
%\references{}
\author{Xuefei Mi}
\note{No further comment}
\seealso{\code{\link{multistageoptimal.grid}}}
\examples{

corr=matrix( c(1,       0.3508,0.3508,0.4979,
               0.3508  ,1,     0.3016,0.5630,
               0.3508,  0.3016,1     ,0.5630,
               0.4979,  0.5630,0.5630,1), 
              nrow=4  
)

multistageoptimal.nlm(N.upper=rep(100,3), corr=corr, Budget=200, CostC=0.5, N.fs=5)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
