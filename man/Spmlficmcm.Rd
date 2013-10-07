\name{Spmlficmcm}
\alias{Spmlficmcm}
\title{
Semiparametric maximum likelihood for interaction in case-mother control-mother
}
\description{
The function builds the nonlinear system from the data, solves the system and assesses the effect of each factor of the model, computes the variance - covariance matrix and deduces from it the standard deviations of each factor.
}
\usage{
Spmlficmcm(fl, N, gmname, gcname, DatfE, typ, vcInp)
}
\arguments{
  \item{fl}{Model formula.
}
  \item{N}{Numeric vector containing eligible number cases and controls in the study population N=(N0, N1).
}
  \item{gmname}{Name of mother genotype variable.
}
  \item{gcname}{Name of offspring genotype variable.
}
\item{DatfE}{\code{data.frame} in long format containing the following variables:\emph{outcome} variable, \emph{mother genotype}, \emph{offspring genotype} and \emph{environmental factors}.
}
\item{typ}{Argument indicating if the data are complete (1) or contain missing offspring genotypes (2).
}
\item{vcInp}{Vector of the initial values of the model parameters.
}
}
\details{
The function \code{Spmlficmcm} builds the nonlinear system from the data and solves the nonlinear system. Then, it uses the log profile likelihood function and the one-step method to estimate the parameters of each factor of the model formula and their standard errors. The Hessian and gradient of the profile likelihood are computed numerically.  The genotype is coded as the number of minor alleles. The model supposes that the distribution of maternal genotype and offspring genotype satisfy the following assumptions: \emph{random mating}, \emph{Hardy-Weinberg equilibrium} and \emph{Mendelian inheritance}. When the data contains missing offspring genotypes, the profile likelihood is summed over the possible genotypes of each child whose genotype is missing. The argument \code{typ} allows the user to specify whether the data is complete or not. Argument \code{vcInp} permits to the user to give the initials values of model parameter. 
Ex: in the following equation log(P/(1-P))=B0+B1*X1+B2*X2+Bm*Gm+Bc*Gc+B2m*X2:Gm, \code{vcIn}=(B0, B1, B2, Bm, Bc, B2m, fp) where fp is the log of the odds of the minor allelic frequency. However, if the user provides no values, the function uses logistic regression to compute the initial B=(B0, B1, B2, Bm, Bc, B2m) and takes 0.1 as the initial value of fp.
}
\value{
 A list containing components
\item{\code{Uim}}{Nonlinear system solution}
\item{\code{MatR}}{Matrix containing the estimates and their standard errors}
\item{\code{Matv}}{Variance - covariance matrix}
}
\references{
Jinbo Chen, Dongyu Lin and Hagit Hochner (2012) Semiparametric Maximum Likelihood Methods for Analyzing Genetic and Environmental Effects with Case-Control Mother-Child Pair Data. \emph{Biometrics} DOI: 10.1111/j.1541-0420.2011.01728.

Moliere Nguile-Makao and Alexandre Bureau (2013), Semi-parametric Maximum likelihood Method for interaction in case-mother control-mother designs: Package SPmlficmcm \emph{under review}.
}

\examples{
# 1-Creation of database
\dontrun{
  set.seed(13200)
  M=20000;
  fl=outc~X1+X2+gm+gnch+X1:gnch+X2:gm;
  theta=0.3
  beta=c(-0.916,0.857,0.588,0.405,-0.693,0.488)
  interc=-2.23
  vpo=c(3,4)
  vprob=c(0.35,0.55)
  vcorr=c(2,1)
  Dataf<-FtSmlrmCMCM(fl,M,theta,beta,interc,vpo,vprob,vcorr)
  rho<-table(Dataf$outc)[2]/20000 # Disease prevalence
         
  # Number of subjects eligible to the study in the population 
  N=c(dim(Dataf[Dataf$outc==0,])[1],dim(Dataf[Dataf$outc==1,])[1])
         
  # Sampling of the study database  
  n0=1232;n1=327; 
  DatfE1<-SeltcEch("outc",n1,n0,"obs",Dataf)


# 2 Creation of missing data on the offspring genotype 
        DatfE=DatfE1 
        gnch<-DatfE["gnch"]
        gnch<-as.vector(gnch[,1])
        gnch1<-sample(c(0,1),length(gnch),replace=TRUE,prob=c(0.91,0.09))
        gnch[gnch1==1]<-NA
        DatfE=DatfE1
        DatfE$gnch<-NULL;DatfE$gnch<-gnch
# 3 Creation of the two databases 
      # DatfEcd :complete data
      # DatfEmd :data with missing genotypes for a subset of children.
        DatfEcd<-DatfE[is.na(DatfE["gnch"])!=TRUE,]
        DatfEmd<-DatfE
        rm(gnch);rm(gnch1) 
# data obtained
DatfEcd[26:30,]
DatfEmd[26:30,]

##4 Estimation of parameters===============================
## model equation         
fl=outc~X1+X2+gm+gnch+X1:gnch+X2:gm;
## Estimation of the parameters (no missing data)
        Rsnm<-Spmlficmcm(fl,N,"gm","gnch",DatfEcd,1)
        #solution of the nonlinear system
        round(Rsnm$Uim,digits=3)
        #estimators
        round(Rsnm$MatR,digits=3)
        #variance covariance matrix
        round(Rsnm$Matv,digits=5)
## Estimation of the parameters (with missing data)
        Rswm<-Spmlficmcm(fl,N,"gm","gnch",DatfEmd,typ=2)
        #solution of the nonlinear system
        round(Rswm$Uim,digits=3)
        #estimates
        round(Rswm$MatR,digits=3)
        #variance covariance matrix
        round(Rswm$Matv,digits=5)
}
}
\keyword{semiparametric}
\keyword{maximum likelihood}
\keyword{case-mother}