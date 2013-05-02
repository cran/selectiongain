\name{multistagecor}
\alias{multistagecor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for calculating correlation matrix in a plant breeding context}
\description{
This function calculates the correlation matrix \eqn{\bm{\Sigma}^{*}} between the test scores X of n stages and the true value y, for n-stage selection in plant breeding. It is n (the number of stages) plus one dimension (for the true value). It is assumed that the genetic structure of the candidates to be selected are genetically fixed, e.g., potential cultivars, clones, inbred lines or testcross progenies of inbred lines with the same or different testers in all stages. The output of this function will be passed to \code{\link{multistagegain}}, which is the main function of this package. The calculation is done for one trait only. For multiple traits, they can be considered as different stages and then calculated by using optimal selection indices (Utz et al. 1994), which is eqivelent to the BLUP estimates of the true value y (Bernardo 1996).

    
}
\usage{
multistagecor(V, L,Rep,index, coe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{V}{Vector of variance components. \eqn{\textbf{V}=\{Vg, Vgl, Vgy, Vgly, Ve \}}, it is refer to the variance among genotypes (\eqn{Vg}), genotype \eqn{\times} location interactions (\eqn{Vgl}), genotype \eqn{\times} year interactions (\eqn{Vgy}), genotype \eqn{\times} year \eqn{\times} location interactions (\eqn{Vgly}) and plot error (\eqn{Ve}). Here, \eqn{Vg} corresponds to \eqn{\sigma_y^2}, which is the genotypic variance of the true genotype, and \eqn{\sigma_y^2} is set equal to 1. Three variances are pre-defined by Longin et al. (2007), they are: "VC1",  \eqn{\textbf{V}=\{1,0.25,0.25,0.5,1\}}; "VC2",  \eqn{\textbf{V}=\{1,0.5,0.5,1,2\}} and "VC3",  \eqn{\textbf{V}=\{1,1,1,2,4\}}.}  
 
  \item{L}{Vector of locations at each stages.}
  
  \item{Rep}{Vector of replicates R at each stage. The values are all equal to one by default.
 }  
 
  \item{index}{Logical. If TRUE (as default), the optimum selection index (Longin et al. 2007), will be used in the calculation. If False, no modification to the test scores X.
 }

  \item{coe}{Logical. If TRUE and index=TRUE, the coefficient of the optimal selection index (BLUP) will be given as an output. In other situations, there will be no output. 
}

}
\details{A special scenario relates to field tests of candidates in plant breeding. Here, the candidates are usually tested in replicated multi-location trials over several years, corresponding to the stages of selection. Thus, besides \eqn{\textbf{N}}, referring to the numbers of candidates to be tested in each stage, the breeder must also decide on the intensity of testing, as reflected by vector for the number of test locations \eqn{\textbf{L}=\{L_1,...,L_n\}} and replications \eqn{\textbf{R}=\{R_1,...,R_n\}}, where \eqn{L_i} and \eqn{R_i} refer to the number of test locations and replications per location, respectively, in stage \eqn{i}. 

If there is no upper limit on \eqn{L_i}, then \eqn{R_i}=1 is optimal for maximizing \eqn{\Delta G} (Longin et al. 2007). Under this scenario, the elements in \eqn{\bm{\Sigma}^{*}} are a rational function of \eqn{\textbf{L}}, \eqn{\textbf{R}} and the vector of variance components \eqn{\textbf{V}}.
 

}
\value{
The output is the correlation matrix \eqn{\bm{\Sigma}^{*}}, which can be used as input parameter of function \code{\link{multistagegain}}.
}
\references{
Bernardo, R. (1996). Best Linear Unbiased Prediction of the Performance of Crosses between Untested Maize Inbreds. \emph{Crop Science}, \bold{4}, 872-876. 

Cochran, W.G. (1951). Improvement by Means of Selection. \emph{Proc. Second Berkeley Symp. on Math. Statist. and Prob.}, 449-470. University of California Press, Berkeley.

Longin, C., Utz, H.F., Reif, J., Wegenast, T., Schipprack, W. and Melchinger, A.E. (2007). Hybrid Maize Breeding with Doubled Haploids: III. Efficiency of Early Testing Prior to Doubled Haploid Production in Two-stage Selection for Testcross Performance. \emph{Theoretical and Applied Genetics}, \bold{115(4)}, 519-527.

Mi, X., Utz, H.F. and Melchinger, A.E. (2007) (submitted). R Package selectiongain: Calculating and Optimizing Expected Gain from Multi-stage Selection \emph{Journal of Statistical Software} 

Utz, H. F., Melchinger, A.E., Seitz, G., Mistele, M. and Zeddies, J. (1994). Economic Aspects of Breeding for Yield and Quality Traits in Forage Maize \emph{Plant Breeding}, \bold{112(2)}, 110-119.

 }
\author{Xuefei Mi}
\note{ No further comments

}
\seealso{ \code{\link{multistagegain}} }
\examples{
Vc2=c(1,0.5,0.5,1,2)
multistagecor(V=Vc2,L=c(2,10))

# "VC2" mentioned in Longin is pre-defined as Vc2=c(1,0.5,0.5,1,2)

multistagecor(V="VC2",L=c(2,10))

# so the above two should have identical results

# further applications
multistagecor(V="VC2",L=c(2,10,20),R=c(1,2,3),index=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Optimization}
% __ONLY ONE__ keyword per line
