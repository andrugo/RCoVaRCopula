CoVaR <- function(y,x,par,par2,dof,gamma,
                  cond.mean,cond.sigma
                  ,dist,type)
{

  alpha<-
  switch (type,
          Gumbel       = Gumbel_CoVaR(level.CoVaR = y , level.VaR=x, par = par),
          RotGumbel    = Rot_Gumbel_CoVaR(level.CoVaR = y , level.VaR = x, par = par),
          Clayton      = Clayton_CoVaR(level.CoVaR = y , level.VaR = x, par = par),
          RotClayton   = Rotated_Clayton_CoVaR(level.CoVaR = y , level.VaR = y, par = par),
          Frank        = Frank_CoVaR(level.CoVaR = y , level.VaR = y, par = par),
          Plackett     = Plackett_CoVaR(level.CoVaR = y , level.VaR = y, par = par),
          Gaussian     = Gaussian_CoVaR(level.CoVaR = y , level.VaR = y, par = par),
          Student      = Student_CoVaR(level.CoVaR = y , level.VaR = y, par = par, par2 = par2),
          BB1          = BB1_CoVaR(level.CoVaR = y , level.VaR = y, par = par, par2 = par2),
          BB7          = BB7_CoVaR(level.CoVaR = y , level.VaR = y, par = par, par2 = par2),
          SCJ          = SCJ_CoVaR(level.CoVaR = y , level.VaR = y, par = par, par2 = par2),
          GumbelUp     = Gumbel_CoVaRUp(level.CoVaR = y , level.VaR=x, par = par),
          RotGumbelUp  = Rot_Gumbel_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par),
          ClaytonUp    = Clayton_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par),
          RotClaytonUp = Rotated_Clayton_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par),
          FrankUp      = Frank_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par),
          PlackettUp   = Plackett_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par),
          GaussianUp   = Gaussian_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par),
          StudentUp    = Student_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par, par2 = par2),
          BB1Up        = BB1_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par, par2 = par2),
          BB7Up        = BB7_CoVaRUp(level.CoVaR = y , level.VaR = y, par = par, par2 = par2),
          SCJUp        = SCJ_CoVaRUp(level.CoVaR = y , level.VaR = x, par = par, par2 = par2)
         )

  perc<-
  switch(dist,
         t     = as.matrix(qt(alpha,df=as.numeric(dof))),
         gauss = as.matrix(qnorm(alpha)),
         tskew = skewtdis_inv(alpha,nu=dof,lambda=gamma)
         )
  CoVaR=NULL
  CoVaR=cond.mean+cond.sigma*as.numeric(perc)
  colnames(CoVaR)<-"CoVaR"
  list(CoVaR=CoVaR,perc=perc,alpha=alpha)
}




