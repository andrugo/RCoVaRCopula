skewtdis_inv <- function(u, nu=nu, lambda=lambda)
{
# PURPOSE: returns the inverse cdf at u of Hansen's (1994) 'skewed t' distribution
#---------------------------------------------------
# USAGE: inv = skewtdis_inv(u,nu,lambda)
# where: U  = a matrix, vector or scalar in the unit interval
#        nu = a matrix or scalar degrees of freedom parameter
#			  lambda = a maxtrix or scalar skewness parameter
#---------------------------------------------------
# RETURNS:
#        a matrix of the inverse cdf at each element of u
# --------------------------------------------------
# SEE ALSO: tdis_cdf, tdis_rnd, tdis_inv, tdis_prb, skewtdis_pdf
#---------------------------------------------------
#
# Based on tdis_inv.m from the "Spatial Econometrics"
# toolbox of James P. LeSage
# http://www.spatial-econometrics.com/
#
#  Andrew Patton
#
#  25 June, 2001

# This code was used in:
#
#	Patton, Andrew J., 2002, "On the Out-of-Sample
# Importance of Skewness and Asymmetric Dependence for
# Asset Allocation", working paper, Department of Economics,
# University of California, San Diego.
#-----------------------------------------------
# Adapted to R code by Andrea Ugolini
#
u=as.matrix(u)
N=NROW(u) # in Matlab = T
k=NCOL(u)

nu=matrix(nu,nrow=N,ncol=1)
lambda=matrix(lambda,nrow=N,ncol=1)

c = gamma((nu+1)/2)/(sqrt(pi*(nu-2))*gamma(nu/2))
a = 4*lambda*c*((nu-2)/(nu-1))
b = sqrt(1 + 3*lambda^2 - a^2)


f1= u < (1-lambda)/2
f2= u >= (1-lambda)/2

inv1 = (1-lambda[f1])/b[f1]*sqrt((nu[f1]-2)/nu[f1])*qt(u[f1]/(1-lambda[f1]),nu[f1])-a[f1]/b[f1]
inv2 = (1+lambda[f2])/b[f2]*sqrt((nu[f2]-2)/nu[f2])*qt(0.5+1/(1+lambda[f2])*(u[f2]-(1-lambda[f2])/2),nu[f2])-a[f2]/b[f2]

inv=matrix(0,nrow=N,ncol=1)
inv[f1]=inv1
inv[f2]=inv2

inv
}
