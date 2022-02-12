\name{multistageoptimum.searchThreeS}
\alias{multistageoptimum.searchThreeS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for optimizing four-stage selection in plant breeding with one marker-assisted selection stage and three phenotypic selection stages}

\description{ This function is used to calculate the maximum of \eqn{\Delta G} based on correlation matrix, which depends on locations, testers and replicates, with a grid search algorithm. The changing correlation matrix of four-stage selection are the testcross progenies of DH lines in one marker-assisted selection (MAS) stage and three phenotypic selection (PS) stages.
}
\usage{
multistageoptimum.searchThreeS (maseff=0.4, VGCAandE,
  VSCA, CostProd, CostTest,  Nf, Budget, N2grid,
  N3grid, N4grid, L2grid, L3grid, L4grid, T2grid,
  T3grid, T4grid, R2, R3, R4, alg,
  detail, fig,alpha.nursery,cost.nursery,
  t2free,parallel.search,saveresult)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{maseff}{ is the efficiency of MAS, if set to NA no marker assited selection or genomic selection is developed in the first stage}


  \item{VGCAandE}{is the vector of variance components of genetic effect, genotype \eqn{\times} location interaction, genotype \eqn{\times} year interaction, genotype \eqn{\times} location \eqn{\times} year interaction and the plot error. When \code{VSCA} is specified, it refers to the general combining ability, otherwise it stands for genetic effect.  The default value is {1,1,1,1,1}. Variances types listed in Longin et al. (2007) can be used. E.g., \code{VGCAandE="VC2"} will set the value as {1,0.5,0.5,1,2}. }

   \item{VSCA}{is the vector of variance components for specific combining ability.}

\item{CostProd}{ contains the initial costs of producing or identifying a candidate in each stage, then the vector should be of lenght four.}

\item{CostTest}{ contains a vector with length n reflecting the cost of evaluating a candidate in the tests performed at stage i, i=1,...,n. The cost might vary in different stages. For this function n=4}

\item{Nf}{ is the number of finally selected candidates.}


\item{Budget}{contains the value of total budget.}


 \item{N2grid}{is the vector of lower and upper limits as well as the grid width of number of candidates in the first field test stage.}

\item{N3grid}{is the vector of lower and upper limits as well as the grid width of number of candidates in the second field test stage.}

\item{N4grid}{is the vector of lower and upper limits as well as the grid width of number of candidates in the third field test stage.}

\item{L2grid}{is the vector of lower and upper limits of number of location as well as the width in the first field test stage.}

\item{L3grid}{is the vector of lower and upper limits of number of location as well as the width in the second field test stage.}

\item{L4grid}{is the vector of lower and upper limits of number of location as well as the width in the third field test stage.}

\item{T2grid}{is the vector of lower and upper limits of number of tester as well as the width in the first field test stage.}

\item{T3grid}{is the vector of lower and upper limits of number of tester as well as the width in the second field test stage.}

\item{T4grid}{is the vector of lower and upper limits of number of tester as well as the width in the third field test stage.}

\item{R2}{is the number of replications in the first field test stage. By default it is 1. }

\item{R3}{is the number of replications in the second field test stage. By default it is 1.}

\item{R4}{is the number of replications in the third field test stage. By default it is 1.}

\item{alg}{is used to switch between two algorithms. If \code{alg = GenzBretz()}, which is by default, the quasi-Monte Carlo algorithm from Genz et al. (2009, 2013), will be used. If \code{alg = Miwa()}, the program will use the Miwa algorithm (Mi et al., 2009), which is an analytical solution of the MVN integral. Miwa's algorithm has higher accuracy (7 digits) than quasi-Monte Carlo algorithm (5 digits). However, its computational speed is slower. We recommend to use the Miwa algorithm.}

\item{detail}{ is the control parameter to decide if the result of all the grids will be given (\code{=TRUE}) or only the maximum (\code{=FALSE}).}


\item{fig}{is the control parameter to decide if a contour plot will be saved in the default folder of R. The default value is \code{FALSE}, which means no figure will be saved.}

\item{alpha.nursery}{a value that should be 0<x<1. The alpha fraction, or amount of genotypes preliminary selected in nurseries, correspond to the fraction entering stage 1 (when MAS is used) or stage 2 (when there is no MAS). It is setted to 1 as default, i.e. no preliminary test "nursery stage". }

\item{cost.nursery}{a vector of length two c([cost of producing a DH line],[cost of testing a DH in nursery]). The default value is {0,0}. }

\item{t2free}{ is a logical value. If =FALSE, the cost of using T4, T3 and T2 testers will be accounted seperately. If =TRUE, the cost of using T4, T3 and T2 testers will be accounted according to number of testers, i.e., CostProd=c(CostProd[1],CostProd[2]*T2,CostProd[3]*(T3-T2),CostProd[4]*(T4-T3) }

\item{parallel.search}{is a logical variable to desided if the multiple cores can be used for computing, by default is FALSE. The users have to notice that assign cores also cost time. So this procedure can only be efficient if the dim >5.}

\item{saveresult}{is a logical variable to save resultfile in saveresult.csv.}

}
\details{

Some breeding programs require more than two phenotypic selection stages. In this programs, a large number of genotypes are assessd for the target trait only in few locations in the first stage and strong selection preasure is applyed. The second and third stages of phenotypic selection are developed in a large number of locations including only a reduced number of genotypes. Even if this stragegy could lead to a reduced selection gain, it could be of major advantage when breeding programs have biological or operative restrictions to conduct large experiments a in large number of locations. This function allows breeders to estimate the possible increase or reduction of selection gain when moving from two stages of phenotypic selection to three stages and also when a rectricted number of genotypes and locations in each of the three stages of phenotypic selection is used.


for the new added to parameters "alpha.nursery" and "cost.nursery" since v2.0.47:

 After producing new DH lines, breeders do NOT go directly for a selection stage in the field, neither for genomic selection. Most of the times, they prefer to make a small field experiment (called "nursery") in which all DH lines are observed and discarded for other traits as disease resistance. That means, all DH lines with poor resistance will be discarded. At the end of the nursery stage only certain amount of DH lines (alpha) advance to the first selection stage (phenotypic or genomic). Specially in maize that makes sense, because in experience around 90 percent of the new DH lines are very weak in terms of per se performance what make them not suitable as new hybrid parents. Then, budget should not be used to make genotyping on or testcrossing with them. Only the alpha fraction should be used for entering the stage 1 of the multistageoptimum.search function.

More details are available in the Crop Science and Computational Statistics papers.

}
\value{
If \eqn{\texttt{detail}} = FALSE, the output of this function is a vector of the optimum allocation i.e., which achieves the maximum \eqn{\Delta G}. Otherwise, the result for all the grid points, which have been calculated, will be exported as a table in the Rgui.}
\references{

A. Genz and F. Bretz. Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195, Springer-Verlag, Heidelberg, 2009.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl and T. Hothorn. mvtnorm: Multivariate normal and t distributions. R package version 0.9-9995, 2013.

E.L. Heffner, A.J. Lorenz, J.L. Jannink, and M.E. Sorrells.  Plant breeding with genomic selection: gain per unit time and cost. Crop Sci. 50: 1681-1690, 2010.

X. Mi, T. Miwa and T. Hothorn. Implement of Miwa's analytical algorithm of multi-normal distribution. R Journal, 1:37-39, 2009.



 }
\author{Jose Marulanda, Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{

VCGCAandError=c(0.4,0.2,0.2,0.4,2)
VCSCA=c(0.2,0.1,0.1,0.2)

#Budget is reduced to 1000 to save computation time

multistageoptimum.searchThreeS(maseff=NA, VGCAandE=VCGCAandError, VSCA=VCSCA,
   alpha.nursery = 0.25, cost.nursery = c(1,0.3), CostProd=c(0,4,4,4), CostTest=c(0,1,1,1),
   Nf=3, Budget=1000, N2grid=c(50,200,50),N3grid=c(10,50,5), N4grid=c(10,20,5),
   L2grid=c(1,2,1), L3grid=c(2,3,1), L4grid=c(4,5,1),
   T2grid=c(1,2,1), T3grid=c(2,3,1), T4grid=c(4,5,1),
   R2=1, R3=1, R4=1, alg=Miwa(), detail=FALSE, fig= FALSE, t2free=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
