Rot_Gumbel_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(level.VaR+x-1+exp(-((-log(1-level.VaR))^par[i]+(-log(1-x))^par[i])^(1/par[i]))))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
  }
  cond
}


Clayton_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(((level.VaR^(-par[i]))+(x^(-par[i]))-1)^(-(1/par[i]))))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
  cond
}

Gumbel_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(exp(-((-log(level.VaR))^par[i]+(-log(x))^par[i])^(1/par[i]))))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
  cond
}

BB1_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par,par2=par2)
{
  cond=NULL
  par=as.matrix(par)
  par2=as.matrix(par2)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+((1+((((level.VaR^(-par[i]))-1)^par2[i])+(((x^(-par[i]))-1)^par2[i]))^(1/par2[i]))^(-1/par[i])))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
  cond
}


Plackett_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(1/(2*(par[i]-1))*(1+(par[i]-1)*(level.VaR+x)-sqrt(((1+(par[i]-1)*(level.VaR+x))^2)-4*par[i]*(par[i]-1)*level.VaR*x))))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
  cond
}


Frank_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(((-1)/par[i])*log(((1-exp(-par[i]))-(1-exp(-par[i]*level.VaR))*(1-exp(-par[i]*x)))/(1-exp(-par[i])))))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
  cond
}

Rotated_Clayton_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(level.VaR+x-1+(((1-level.VaR)^-par[i])+((1-x)^-par[i])-1)^(-1/par[i])))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
 cond
}


BB7_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par,par2=par2)
{
  cond=NULL
  par=as.matrix(par)
  par2=as.matrix(par2)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) (1-level.VaR-x+(1-(1-((1-(1-level.VaR)^par[i])^(-par2[i])+(1-(1-x)^par[i])^(-par2[i])-1)^(-1/par2[i]))^(1/par[i])))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
}
  cond
}


Student_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par,par2=par2)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) pCopula(c(x,1-level.VaR), tCopula(as.numeric(par[i]),dim=2,df=as.integer(par2)))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,0, tol = 10^-16)$x
}
  out=1-cond
}

Gaussian_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par)
{
  cond=NULL
  par=as.matrix(par)
  N=length(par)
  for (i in 1:N){
    CoVaR=function(x) pCopula(c(x,1-level.VaR), normalCopula(par[i]))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,0, tol = 10^-16)$x
}
  out=1-cond
}

SCJ_CoVaRUp<-function(level.CoVaR=level.CoVaR,level.VaR=level.VaR,par=par,par2=par2)
{
  #par = TauUpper
  #par2 = TauLower
  cond=NULL
  par=as.matrix(par)
  par2=as.matrix(par2)
  N=length(par)
  for (i in 1:N){
    K =  1/log2(2-par[i]);
    G = -1/log2(par2[i]);
    K1 =  1/log2(2-par2[i]); # switching the upper and lower measures
    G1 = -1/log2(par[i]);
    CoVaR=function(x) (1-level.VaR-x+(0.5*((1-((1-(((1-((1-level.VaR)^K))^(-G))+((1-((1-x)^K))^(-G))-1)^(-1/G))^(1/K)))
                                           +(1-(1-level.VaR)) + (1-(1-x)) - 1 + 1-((1-(((1-((1-(1-level.VaR))^K1))^(-G1))
                          +((1-((1-(1-x))^K1))^(-G1))-1)^(-1/G1))^(1/K1)))))-((1-level.CoVaR)*(1-level.VaR))
    cond[i]=fzero(CoVaR,1, tol = 10^-16)$x
  }
  cond
}

