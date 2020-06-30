\name{MCMCglmm}
\alias{MCMCglmm}
\title{Multivariate Generalised Linear Mixed Models}
\description{
 Markov chain Monte Carlo Sampler for Multivariate Generalised Linear Mixed 
Models with special emphasis on correlated random effects arising from pedigrees 
and phylogenies (Hadfield 2010). Please read the course notes: \code{vignette("CourseNotes", 
"MCMCglmm")} or the overview \code{vignette("Overview", "MCMCglmm")}
}
\usage{
MCMCglmm(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, 
    data,start=NULL, prior=NULL, tune=NULL, pedigree=NULL, nodes="ALL",
    scale=TRUE, nitt=13000, thin=10, burnin=3000, pr=FALSE,
    pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=TRUE,
    saveZ=TRUE, saveXL=TRUE, slice=FALSE, ginverse=NULL, trunc=FALSE)
}
\arguments{
  \item{fixed}{\code{\link{formula}} for the fixed effects, multiple responses 
are passed as a matrix using cbind}
  \item{random}{\code{\link{formula}} for the random effects. Multiple random terms can be passed using the \code{+} operator, and in the most general case each random term has the form \code{variance.function(formula):linking.function(random.terms)}. Currently, the only \code{variance.functions} available are \code{idv}, \code{idh}, \code{us}, \code{cor[]} and \code{ante[]}. \code{idv} fits a constant variance across all components in \code{formula}. Both \code{idh} and \code{us} fit different variances across each component in \code{formula}, but \code{us} will also fit the covariances. \code{corg} fixes the variances along the diagonal to one and \code{corgh} fixes the variances along the diagonal to those specified in the prior.  \code{cors} allows correlation submatrices. 
\code{ante[]} fits ante-dependence structures of different order (e.g ante1, ante2), and the number can be prefixed by a \code{c} to hold all regression coefficients of the same order equal. The number can also be suffixed by a \code{v} to hold all innovation variances equal (e.g \code{antec2v} has 3 parameters).  The \code{formula} can contain both factors and numeric terms (i.e. random regression) although it should be noted that the intercept term is suppressed.  The (co)variances are the (co)variances of the \code{random.terms} effects.  Currently, the only \code{linking.functions} available are \code{mm} and \code{str}. \code{mm} fits a multimembership model where multiple random terms are separated by the \code{+} operator. \code{str} allows covariances to exist between multiple random terms that are also separated by the \code{+} operator. In both cases the levels of all multiple random terms have to be the same.  For simpler models the \code{variance.function(formula)} and \code{linking.function(random.terms)} can be omitted and the model syntax has the simpler form \code{~random1+random2+...}.  There are two reserved variables: \code{units} which index rows of the response variable and \code{trait} which index columns of the response variable} 
  \item{rcov}{\code{\link{formula}} for residual covariance structure. This has to be set up so that each data point is associated with a unique residual. For example a multi-response model might have the R-structure defined by \code{~us(trait):units}}
  \item{family}{optional character vector of trait distributions. Currently, 
\code{"gaussian"}, \code{"poisson"},  \code{"categorical"},  
\code{"multinomial"},  \code{"ordinal"}, \code{"threshold"}, \code{"exponential"}, \code{"geometric"}, \code{"cengaussian"},  
\code{"cenpoisson"},  \code{"cenexponential"},  \code{"zipoisson"}, \code{"zapoisson"}, \code{"ztpoisson"}, \code{"hupoisson"}, \code{"zibinomial"}, \code{"threshold"} and \code{nzbinom} are supported, where the prefix \code{"cen"} means censored, the prefix \code{"zi"} means zero inflated, the prefix \code{"za"} means zero altered, the prefix \code{"zt"} means zero truncated and the prefix \code{"hu"} means hurdle. If \code{NULL}, \code{data} needs to contain a 
\code{family} column.}
  \item{mev}{optional vector of measurement error variances for each data point 
for random effect meta-analysis.}
  \item{data}{\code{data.frame}}
  \item{start}{optional list having 4 possible elements: 
\code{R} (R-structure) \code{G} (G-structure) and \code{liab} (latent variables or liabilities) should contain the starting values where \code{G} itself is also a list with as many elements as random effect components. The fourth element \code{QUASI} should be logical: if \code{TRUE} starting latent variables are obtained heuristically, if \code{FALSE} then they are sampled from a Z-distribution}
  \item{prior}{optional list of prior specifications having 3 possible elements: 
\code{R} (R-structure) \code{G} (G-structure) and \code{B} (fixed effects). \code{B} is a list containing the expected value (\code{mu}) and a 
(co)variance matrix (\code{V}) representing the strength of belief: the defaults are \code{B$mu}=0 and \code{B$V}=I*1e+10, where where I is an identity matrix of appropriate dimension.  The priors for the variance structures (\code{R} and \code{G}) are lists with the expected (co)variances (\code{V}) and degree of belief parameter (\code{nu}) for the inverse-Wishart, and also the mean vector (\code{alpha.mu}) and covariance matrix (\code{alpha.V}) for the redundant working parameters. The defaults are \code{nu}=0, \code{V}=1, \code{alpha.mu}=0, and \code{alpha.V}=0. When \code{alpha.V} is non-zero, parameter expanded algorithms are used.} 
  \item{tune}{optional (co)variance matrix defining the proposal distribution 
for the latent variables. If NULL an adaptive algorithm is used which ceases to 
adapt once the burn-in phase has finished.}
  \item{pedigree}{ordered pedigree with 3 columns id, dam and sire or a 
\code{phylo} object. This argument is retained for back compatibility - see ginverse argument for a more general formulation.}
  \item{nodes}{pedigree/phylogeny nodes to be estimated. The default, 
\code{"ALL"} estimates effects for all individuals in a pedigree or nodes in a 
phylogeny (including ancestral nodes). For phylogenies \code{"TIPS"} estimates 
effects for the tips only, and for pedigrees a vector of ids can be passed to 
\code{nodes} specifying the subset of individuals for which animal effects are 
estimated.  Note that all analyses are equivalent if omitted nodes have missing 
data but by absorbing these nodes the chain max mix better. However, the 
algorithm may be less numerically stable and may iterate slower, especially for 
large phylogenies.}
  \item{scale}{logical: should the phylogeny (needs to be ultrametric) be scaled 
to unit length (distance from root to tip)?}
  \item{nitt}{number of MCMC iterations}
  \item{thin}{thinning interval}
  \item{burnin}{burnin}
  \item{pr}{logical: should the posterior distribution of random effects be 
saved?}
  \item{pl}{logical: should the posterior distribution of latent variables be 
saved?}
  \item{verbose}{logical: if \code{TRUE} MH diagnostics are printed to screen}
  \item{DIC}{logical: if \code{TRUE} deviance and deviance information criterion are calculated}
  \item{singular.ok}{logical: if \code{FALSE} linear dependencies in the fixed effects are removed. if \code{TRUE} they are left in an estimated, although all information comes form the prior}
  \item{saveX}{logical: save fixed effect design matrix}
  \item{saveZ}{logical: save random effect design matrix}
  \item{saveXL}{logical: save structural parameter design matrix}
  \item{slice}{logical: should slice sampling be used? Only applicable for binary traits with independent residuals}
  \item{ginverse}{a list of sparse inverse matrices (\eqn{{\bf A^{-1}}}{solve(A)}) that are proportional to the covariance structure of the random effects. The names of the matrices should correspond to columns in \code{data} that are associated with the random term. All levels of the random term should appear as rownames for the matrices.}
    \item{trunc}{logical: should latent variables in binary models be truncated to prevent under/overflow (+/-20 for categorical/multinomial models and +/-7 for threshold/probit models)? }
}
\value{
  \item{Sol}{Posterior Distribution of MME solutions, including fixed effects}
  \item{VCV}{Posterior Distribution of (co)variance matrices}
  \item{CP}{Posterior Distribution of cut-points from an ordinal model}
  \item{Liab}{Posterior Distribution of latent variables}
  \item{Fixed}{list: fixed formula and number of fixed effects}
  \item{Random}{list: random formula, dimensions of each covariance matrix, number of levels per covariance matrix, and term in random formula to which each covariance belongs}
  \item{Residual}{list: residual formula, dimensions of each covariance matrix, number of levels per covariance matrix, and term in residual formula to which each covariance belongs}
  \item{Deviance}{deviance -2*log(p(y|...))}
  \item{DIC}{deviance information criterion}
  \item{X}{sparse fixed effect design matrix}
  \item{Z}{sparse random effect design matrix}
  \item{XL}{sparse structural parameter design matrix}
  \item{error.term}{residual term for each datum}
  \item{family}{distribution of each datum}
  \item{Tune}{(co)variance matrix of the proposal distribution for the latent variables}
  \item{meta}{logical; was \code{mev} passed?}
}
\references{
	General analyses: Hadfield, J.D. (2010) Journal of Statistical Software 33 2 1-22
	
	Phylogenetic analyses: Hadfield, J.D. & Nakagawa, S. (2010) Journal of Evolutionary Biology 23 494-508
	
	Background Sorensen, D. & Gianola, D. (2002) Springer
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link[coda]{mcmc}}}
\examples{

# Example 1: univariate Gaussian model with standard random effect
 
data(PlodiaPO)  
model1<-MCMCglmm(PO~1, random=~FSfamily, data=PlodiaPO, verbose=FALSE,
 nitt=1300, burnin=300, thin=1)
summary(model1)

# Example 2: univariate Gaussian model with phylogenetically correlated
# random effect

data(bird.families) 

phylo.effect<-rbv(bird.families, 1, nodes="TIPS") 
phenotype<-phylo.effect+rnorm(dim(phylo.effect)[1], 0, 1)  

# simulate phylogenetic and residual effects with unit variance

test.data<-data.frame(phenotype=phenotype, taxon=row.names(phenotype))

Ainv<-inverseA(bird.families)$Ainv

# inverse matrix of shared phyloegnetic history

prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

model2<-MCMCglmm(phenotype~1, random=~taxon, ginverse=list(taxon=Ainv),
 data=test.data, prior=prior, verbose=FALSE, nitt=1300, burnin=300, thin=1)

plot(model2$VCV)

}
\keyword{models}


