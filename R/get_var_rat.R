 get_var_rat<-function(var.in1,var.in2,var.dof1,var.dof2){
 # get unbiased estimates of the mean variance ratio by accounting for the degrees of freedom of the denominator
#  https://people.richland.edu/james/lecture/m113/f_test.html
varratio<-var.in1/var.in2*(var.dof2-2)/var.dof2
return(varratio)
}
 