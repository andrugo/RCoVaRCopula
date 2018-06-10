# CoVaR value evaluting at quantile CoVaR and quantile VaR.

# description
Calculate the conditional quantile or CoVaR with different type of Copula and marginal distribution. In this package several bivariate copula families are included for bivariate analysis.
It provides functionality of elliptical (Gaussian and Student t) as well as Archimedean (Clayton, Gumbel, Frank, Plackett, BB1, SCJ, rotated clayton and rotated Gumbel) copulas to cover a large bandwidth of possible dependence structures.



# Author 
Andrea Ugolini <andreaugolini@me.com> \\
Juan Carlos Reboredo Noguiera <juancarlos.reboredo@usc.es>

# References
Reboredo, J. C., & Ugolini, A. (2016). 
Quantile dependence of oil price movements and stock returns. Energy Economics, 54, 33-49.

# Example
#### RCoVaRCopula
load("Data_demo.Rdata")
source("CoVaR.R")
source("DynCopulaCoVaR.R")
source("DynCopulaCoVaRUpper.R")
source("skewtdis_inv.R")
require("pracma")
require("copula")

#### CoVaR Downside

CoVaR1part = CoVaR(0.05,0.05,par=par1_1,par2=par2_1,dof=tailBrazil,gamma=asyBrazil,cond.mean=meanBrasil1,
             cond.sigma=sigmaBrasil1,dist="tskew",type="Student")

CoVaR2part = CoVaR(0.05,0.05,par=par1_2,par2=par2_2,dof=tailBrazil,gamma=asyBrazil,cond.mean=meanBrasil2,
             cond.sigma=sigmaBrasil2,dist="tskew",type="Student")
#### CoVaR Upside

CoVaR1partUp = CoVaR(0.95,0.95,par=par1_1,par2=par2_1,dof=tailBrazil,gamma=asyBrazil,cond.mean=meanBrasil1,
               cond.sigma=sigmaBrasil1,dist="tskew",type="StudentUp")

CoVaR2partUp = CoVaR(0.95,0.95,par=par1_2,par2=par2_2,dof=tailBrazil,gamma=asyBrazil,cond.mean=meanBrasil2,
               cond.sigma=sigmaBrasil2,dist="tskew",type="StudentUp")

CoVaR1D = CoVaR1part$CoVaR
CoVaR2D = CoVaR2part$CoVaR
CoVaR1U = CoVaR1partUp$CoVaR
CoVaR2U = CoVaR2partUp$CoVaR

TimeCoVaRD = rbind(CoVaR1D,CoVaR2D)
TimeCoVaRU = rbind(CoVaR1U,CoVaR2U)

#### Plot
plot(as.matrix(TimeCoVaRD),type="l",col="blue",
     ylim=c(-0.5,0.5),xlab="Time",ylab="")
lines(VaR,col="black",lty=2)
lines(TimeCoVaRU,col="red",lty=4)
lines(VaRup,col="green",lty=3)
abline(h=0,col="gray33")

R Code CoVaR with Copula
