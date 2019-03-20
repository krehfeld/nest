#' Confidence intervals for variances, and ratios of variances.
#' 
#' Confidence intervals for variances, and ratios of variances
#' 
#' 
#' @aliases var_ci varrat_ci get_var_rat
#' @usage var_ci(var.in, var.dof, pval = 0.1) varrat_ci(varratio, var.dof1,
#' var.dof2, pval = 0.1) get_var_rat(var.in1,var.in2,var.dof1,var.dof2)
#' @param var.in Variance estimate (\code{var_ci})
#' @param var.dof Degrees of freedom (\code{var_ci})
#' @param varratio Ratio of the two variances (\code{varrat_ci})
#' @param var.dof1 Degrees of freedom for variance ratio (nominator)
#' @param var.dof2 Degrees of freedom for variance ratio (denominator)
#' @param pval P-value for which the CIs are desired
#' @param var.in1 Nominator variance
#' @param var.in2 Denominator variance
#' @return Confidence intervals for variances and variance ratios.
#' @author Kira Rehfeld and Thomas Laepple
#' @seealso \code{\link{tsc_dep_var}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' tmp<-var_ci(1,10,pval=0.1)
#' varratio<-get_var_rat(var.in1=0.9,var.in2=1,var.dof1=10,var.dof2=10)
#' tmp2<-varrat_ci(varratio,10,10,pval=0.1)
#' 
#' 
#' @export var_ci
var_ci<-function(var.in,var.dof,pval=0.1)
{

 lim.1<-var.in*qchisq(c(1-pval/2),var.dof)/(var.dof)
 lim.2<-var.in*qchisq(c(pval/2),var.dof)/(var.dof)
res<-list(lo=lim.2,up=lim.1)
 return(res)
}

