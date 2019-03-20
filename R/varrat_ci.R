varrat_ci<-function(varratio,var.dof1,var.dof2,pval=0.1)
#Confidence intervals for Ratios based on the F-distribution
#Input: varratio, and the dofs...
#Output: lower and upper confidence intervals
{ 
if (!(length(varratio)==length(var.dof1))&&(!length(varratio)==length(var.dof2))) {stop("same lengths must be provided")}
if (!(is.numeric(varratio))||!(is.numeric(var.dof1))||!is.numeric(var.dof2)) {stop("non-numeric arguments")}

res<-matrix(NA,nrow=length(varratio),ncol=2)

for (i in 1:length(varratio)){
# print(i)
# print(var.dof1[i])
# print(var.dof2[i])
# print(varratio[i])
QF<-qf(p=c(pval/2,(1-pval/2)),df1=var.dof1[i],df2=var.dof2[i])
tmp<-QF*varratio[i]
#print(tmp)    
res[i,]<-tmp

     }
     return(res)
}
