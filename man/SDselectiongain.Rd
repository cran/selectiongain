\name{SDselectiongain}
\alias{SDselectiongain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function for calculating the standrd deviation of selection gain
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to calculate the standard deviation of sel gain acording to longin 2015
}
\usage{
SDselectiongain(Ob, maseff, VGCAandE, VSCA, VLine, years, Genotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ob}{
%%     ~~Describe \code{Ob} here~~
matrix object produced by the function multistageoptimum.search or multistageoptiumum.grid}
  \item{maseff}{
%%     ~~Describe \code{maseff} here~~
is the efficiency of marker-assisted selection (MAS). The default value is NA, which means there is no MAS. If a value between 0 and 1 is assigned to \code{maseff}, then the first selection stage will be considered as MAS (Heffner et al., 2010). The value of MAS is recommanded to be higher than 0.1 to avoid illshaped correlation matrix.}
  \item{VGCAandE}{
%%     ~~Describe \code{VGCAandE} here~~
is the vector of variance components of genetic effect, genotype \eqn{\times} location interaction, genotype \eqn{\times} year interaction, genotype \eqn{\times} location \eqn{\times} year interaction and the plot error. When \code{VSCA} is specified, the VGCAandE refers to the general combining ability, otherwise it stands for genetic effect.  The default value is {1,1,1,1,1}. Variances types listed in Longin et al. (2007) can be used. For example, \code{VGCAandE="VC2"} will set the value as {1,0.5,0.5,1,2}.}
  \item{VSCA}{
%%     ~~Describe \code{VSCA} here~~
is the vector of variance components for specific combining ability. The default value is {0,0,0,0}. }
  \item{VLine}{
%%     ~~Describe \code{VLine} here~~
is the vector of variance components for line per se. The default value is {0,0,0,0,0}.}
  \item{years}{
%%     ~~Describe \code{years} here~~
Duration of the breeding scheme in years, it is used only to compute the anual selection gain }
  \item{Genotypes}{
%%     ~~Describe \code{Genotypes} here~~
character vector to indicate the function which variance components we are using. Pssible values are "Hybrids" if we are using GCA and SCA variance components or "Lines" if we are using line perse variance components}}
  \details{
%%  ~~ If necessary, more details than the description above ~~

for the new added to parameters "alpha.nursery" and "cost.nursery" since v2.0.47:

 After producing new DH lines, breeders do NOT go directly for a selection stage in the field, neither for genomic selection. Most of the times, they prefer to make a small field experiment (called "nursery") in which all DH lines are observed and discarded for other traits as disease resistance. That means, all DH lines with poor resistance will be discarded. At the end of the nursery stage only certain amount of DH lines (alpha) advance to the first selection stage (phenotypic or genomic). Specially in maize that makes sense, because in experience around 90 percent of the new DH lines are very weak in terms of per se performance what make them not suitable as new hybrid parents. Then, budget should not be used to make genotyping on or testcrossing with them. Only the alpha fraction should be used for entering the stage 1 of the multistageoptimum.search function.

More details are available in the Crop Science and Computational Statistics papers.}
  \value{
%%  ~Describe the value returned
The output is equivalent to the matrix object produced by the functions multistageoptimum.search or multistageoptimum.grid but with two columns added, one for the values of the anual selection gain and the second for the standard deviation of selection gain}
  \references{
C. Longin, X. Mi and T. Wuerschum. Genomic selection in wheat: optimum allocation of test resources and comparison of breeding strategies for line and hybrid breeding. Theoretical and Applied Genetics 128: 1297-1306. 2015.

C. Longin, H.F. Utz., J. Reif, T. Wegenast, W. Schipprack and A.E. Melchinger. Hybrid maize breeding with doubled haploids: III. Efficiency of early testing prior to doubled haploid production in two-stage selection for testcross performance. Theor. Appl. Genet. 115: 519-527, 2007.

E.L. Heffner, A.J. Lorenz, J.L. Jannink, and M.E. Sorrells.  Plant breeding with genomic selection: gain per unit time and cost. Crop Sci. 50: 1681-1690, 2010.}

  \author{
%%  ~~who you are~~
Jose Marulanda}
  \note{
%%  ~~further notes~~
no further comment}

%% ~Make other sections like Warning with \section{Warning }{....} ~

  \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
selectiongain()}
  \examples{

CostProd =c(0.5,1,1)
CostTest = c(0.5,1,1)
Budget=1021
# Budget is very small here to save time in package checking
# for the example in Heffner's paper, please change it to Budget=10021

VCGCAandError=c(0.4,0.2,0.2,0.4,2)
VCSCA=c(0.2,0.1,0.1,0.2)
Nf=10
maseff=0.4
years=7
# this breeding scheme takes 7 years from the initial cross to the final field testing.
# See references for more details


Ob<-multistageoptimum.search (maseff=maseff, VGCAandE=VCGCAandError,
VSCA=VCSCA, CostProd = CostProd, CostTest = CostTest,
Nf = Nf, Budget = Budget, N2grid = c(11, 1211, 30),
N3grid = c(11, 211, 5), L2grid=c(1,1,1), L3grid=c(6,6,1),
T2grid=c(1,2,1), T3grid=c(3,5,1), R2=1, R3=1, alg = Miwa(),
detail=TRUE, fig=FALSE, t2free=TRUE)

SDselectiongain(Ob=Ob,maseff=maseff,VGCAandE=VCGCAandError,VSCA=VCSCA,
                years=years,Genotypes="Hybrids")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ Optimization }
% __ONLY ONE__ keyword per line
