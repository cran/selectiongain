\name{multistagecor}
\alias{multistagecor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating correlation matrix in a plant breeding context}
\description{
This function is used to calculate the (n+1)-dimensional correlation matrix \eqn{\bm{\Sigma}^{*}} of y and X, where y is the true value (genotypic value in plant breeding) and \eqn{\bm{X}=\{X_1,...X_n \}} are the values of y's observations or selection indices, which are linear combinations of the values of observation from each selection stage.

In a plant breeding context, it is assumed that the genetic structure of the candidates to be selected are genetically fixed, e.g., potential cultivars, clones, inbred lines or testcross progenies of inbred lines with the same or different testers in all stages.


}
\usage{
multistagecor(maseff,VGCAandE,VSCA,VLine,ecoweight,rhop,
              T,L,M,Rep,index, indexTrait, covtype, detail,
              VGCAandE2, VSCA2, COVgca, COVsca, maseff2, q12, q22)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{maseff}{ is the efficiency of marker-assisted selection (MAS). The default value is NA, which means there is no MAS. If a value between 0 and 1 is assigned to \code{maseff}, then the first selection stage will be considered as MAS (Heffner et al., 2010). The value of MAS is recommanded to be higher than 0.1 to avoid illshaped correlation matrix.}

  \item{VGCAandE}{is the vector of variance components of genetic effect, genotype \eqn{\times} location interaction, genotype \eqn{\times} year interaction, genotype \eqn{\times} location \eqn{\times} year interaction and the plot error. When \code{VSCA} is specified, the VGCAandE refers to the general combining ability (hybrid breeding), otherwise it stands for genetic effect (line breeding).  The default value is {1,1,1,1,1}. Variances types listed in Longin et al. (2007) can be used. For example, \code{VGCAandE="VC2"} will set the value as {1,0.5,0.5,1,2}. }

  \item{VSCA}{is the vector of variance components for specific combining ability (hybrid breeding). The default value is {0,0,0,0}. }

  \item{VLine}{Only to be used if parental and testcross selection are performed in a breeding strategy, For an example see the paper "Wegenast, Longin... 2008. Hybrid maize breeding with doubled haploids. IV". If this strategy is implemented, then Vline correspond to the vector of variance components for the parents (line per se). The default value is {0,0,0,0,0}.}

  \item{ecoweight}{is the vector of economic weight. In the case of simultaneos selection of two traits, this vector contains two elements, each corresponding to economical weigth of each trait}

  \item{rhop}{is the genetic correlation between line per se performance and GCA}

  \item{T}{is the vector of number of testers at each stage. If there is no tester applied in a certain stage, the value at this stage has to be 1. }

  \item{L}{is the vector of number of locations at each stage.}

  \item{M}{is the vector of tester type, i.e., number of unrelated inbred lines combined in a single tester in stage j.}

  \item{Rep}{is the vector of number of replications at each stage.}

  \item{index}{is the control parameter. If it equals TRUE, the optimum selection index of Longin et al. (2007) will be used in the calculation of correlation matrix without MAS.}

  \item{indexTrait}{is the control parameter for the simultaneous selection of two traits. Possible options are: "Optimum"(default), "Base" and "Restricted" for the implementation of the well known optimum, base and restricted selection indexes in plant breeding.}

   \item{covtype}{is the type of the covariance. Longin's type (\code{covtype}=c("LonginII")) is used by default. For the simultaneous selection of two traits possible covtypes are "2traits_PS", "2traits_GS" , "2traits_GS-PS", "2traits_PS-PS", "2traits_GS-PS-PS". If any of these five option is selected the calculation of correlation matrix will use the variance components of the two traits. If the user also require marker assited selection, the prediction accuracy of MAS for both traits should be also given to the function. Finally, if two traits are selected simultaneously, the desired index have to be defined in indexTrait}

   \item{detail}{ is the control parameter to decide if the correlation matrix, optimal selection index and covariance matrix will be returned (\code{=TRUE}) or only the correlation matrix (\code{FALSE}). The default value is \code{FALSE}.}

   \item{VGCAandE2}{In the case of simultaneos selection of two traits (index selection) it is the vector of variance components of genetic effect, genotype \eqn{\times} location interaction, genotype \eqn{\times} year interaction, genotype \eqn{\times} location \eqn{\times} year interaction and the plot error for the second trait. When \code{VSCA2} is specified, the VGCAandE refers to the general combining ability, otherwise it stands for genetic effect of the second trait.  The default value is {0,0,0,0,0}, meaning no simultaneos selection of two traits.}

  \item{VSCA2}{In the case of simultaneos selection of two traits (index selection) it is the vector of variance components for specific combining ability for the second trait. The default value is {0,0,0,0}. The default value is {0,0,0,0} }

  \item{COVgca}{In the case of simultaneos selection of two traits (index selection) is the vector of covariance components of: genetic effect, genotype \eqn{\times} location interaction, genotype \eqn{\times} year interaction, genotype \eqn{\times} location \eqn{\times} year interaction and the plot error.}

  \item{COVsca}{In the case of simultaneos selection of two traits (index selection) is the vector of covariance components of the specific combining ability effects as follows : sca, sca \eqn{\times} location interaction, sca \eqn{\times} year interaction, sca \eqn{\times} location \eqn{\times} year interaction. .}


   \item{maseff2}{ is the efficiency of marker-assisted selection (MAS) for the second trait. The default value is NA, which means there is no MAS and there is not simultaneous selection of two traits. If a value between 0 and 1 is assigned to \code{maseff2}, then it is assumed that the breeder want to optimize breeding strategies for the simultaneos selection of two traits and also including marker assited selection. In this case, appropiate options have to be selected in covtype and indexTrait. The value of MAS is recommended to be higher than 0.1 to avoid illshaped correlation matrix.}

  \item{q12}{ is the proportion of genetic variance associated with markers for trait 1 as defined by "Dekkers, JCM. 2007. Prediction of response to marker-assited..."" This parameter is only needed in the case of simultaneos selection of two traits (index selection)}

  \item{q22}{ is the proportion of genetic variance associated with markers for trait 2 as defined by "Dekkers, JCM. 2007. Prediction of response to marker-assited..."" This parameter is only needed in the case of simultaneos selection of two traits (index selection)}

}
\value{
The default output is a matrix with dimension n+1 and can be used as input parameter of function multistagegain. When  value of detail=TRUE, the correlation matrix, optimal selection index and covariance matrix will be given. If covtype are set to: "2traits_PS" , "2traits_GS" , "2traits_GS-PS" , "2traits_PS-PS" , or "2traits_GS-PS-PS" , the output will be a list of seven matrices as follows: (1) correlation matrix for the index, (2) estimates of the relative index weights B (betas) for each trait in each stage, (3) covariance matrix for the index (4) correlation matrix for trait 1, (5) correlation matrix for trair 2, (6) matrix of genotypic covariances and (7) matrix of phenotypic covariances
}

\references{

C. Longin, H.F. Utz., J. Reif, T. Wegenast, W. Schipprack and A.E. Melchinger. Hybrid maize breeding with doubled haploids: III. Efficiency of early testing prior to doubled haploid production in two-stage selection for testcross performance. Theor. Appl. Genet. 115: 519-527, 2007.

E.L. Heffner, A.J. Lorenz, J.L. Jannink, and M.E. Sorrells.  Plant breeding with genomic selection: gain per unit time and cost. Crop Sci. 50: 1681-1690, 2010.


 }
\author{ Xuefei Mi}
\note{ no further comment

}
\seealso{ selectiongain() }
\examples{
# example for calculating correlation matrix without MAS
multistagecor(VGCAandE=c(1,0.5,0.5,1,2),L=c(2,10),T=c(1,1),Rep=c(1,1))
multistagecor(VGCAandE="VC2",L=c(2,10),T=c(1,1),Rep=c(1,1),index=TRUE)

# example for calculating correlation matrix with MAS in the first stage
VCgca=c(0.40,0.20,0.20,0.40,2.00)
VCsca=c(0.20,0.10,0.10,0.20)
corr.matrix = multistagecor (maseff=0.40, VGCAandE=VCgca,
VSCA=VCsca, T=c(1,1,5), L=c(1,3,8), Rep=c(1,1,1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
