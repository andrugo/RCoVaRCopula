# RCoVaRCopula
data(Data_demo)

# CoVaR Downside

CoVaR1part = CoVaR(0.05,0.05,par=par1_1,par2=par2_1,dof=tailBrazil,gamma=asyBrazil,cond.mean=meanBrasil1,
             cond.sigma=sigmaBrasil1,dist="tskew",type="Student")

CoVaR2part = CoVaR(0.05,0.05,par=par1_2,par2=par2_2,dof=tailBrazil,gamma=asyBrazil,cond.mean=meanBrasil2,
             cond.sigma=sigmaBrasil2,dist="tskew",type="Student")
# CoVaR Upside

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

# Plot
plot(as.matrix(TimeCoVaRD),type="l",col="blue",
     ylim=c(-0.5,0.5),xlab="Time",ylab="")
lines(VaR,col="black",lty=2)
lines(TimeCoVaRU,col="red",lty=4)
lines(VaRup,col="green",lty=3)
abline(h=0,col="gray33")

R Code CoVaR with Copula
