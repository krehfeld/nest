\name{var_ci}
\alias{var_ci}
\alias{varrat_ci}
\alias{get_var_rat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence intervals for variances, and ratios of variances.}
\description{Confidence intervals for variances, and ratios of variances}
\usage{
var_ci(var.in, var.dof, pval = 0.1)
varrat_ci(varratio, var.dof1, var.dof2, pval = 0.1)
get_var_rat(var.in1,var.in2,var.dof1,var.dof2)
}

\arguments{
  \item{var.in}{Variance estimate (\code{var_ci})}
  \item{var.dof}{Degrees of freedom (\code{var_ci})}
  \item{varratio}{Ratio of the two variances (\code{varrat_ci})}
  \item{var.dof1}{Degrees of freedom for variance ratio (nominator)}
  \item{var.dof2}{Degrees of freedom for variance ratio (denominator)}
  \item{pval}{P-value for which the CIs are desired}
  \item{var.in1}{Nominator variance}
  \item{var.in2}{Denominator variance}


}

\value{Confidence intervals for variances and variance ratios.}

\author{Kira Rehfeld and Thomas Laepple}

\seealso{
\code{\link{tsc_dep_var}}
}
\examples{

tmp<-var_ci(1,10,pval=0.1)
varratio<-get_var_rat(var.in1=0.9,var.in2=1,var.dof1=10,var.dof2=10)
tmp2<-varrat_ci(varratio,10,10,pval=0.1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
