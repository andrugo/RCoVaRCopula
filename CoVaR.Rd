% Please edit documentation in R/CoVaR.R
\name{CoVaR}
\alias{CoVaR}
\title{CoVaR with Copula}
\usage{
CoVaR(y, x, par, par2, dof, gamma, cond.mean, cond.sigma, dist, type)
}
\arguments{
\item{y}{CoVaR confidence livel.}

\item{x}{VaR confidence livel.}

\item{par}{Copula paramenter.}

\item{par2}{Second paramenter for bivariate copulas with two parameters (t, BB1, SCJ).}

\item{dof}{Degree of Freedom if a inverse distribution is a Student-t or Skew-t.}

\item{gamma}{Asymmetry parameter for Skew-t.}

\item{cond.mean}{Conditional mean of y.}

\item{cond.sigma}{Conditional mean of y.}

\item{dist}{Type of inverse distribution function: \cr
      \code{gauss} =  Gaussian distribution function \cr
      \code{t} = Student-t distribution function\cr
      \code{tskew} = Skew-t distribution function\cr
  }    
\item{type}{Type of Copula function: \cr
	   \code{Gaussian} = Gaussian copula \cr
	   \code{Student} = Student-t copula (t-copula) \cr
     \code{Clayton} = Clayton copula \cr
     \code{Gumbel} = Gumbel copula \cr
	   \code{Frank} = Frank copula \cr
	   \code{Plackett} = Plackett copula \cr 
		 \code{BB1} = BB1 copula \cr
		 \code{SCJ} = Symetric Joe-Clayton copula \cr
		 \code{RotClayton} = rotated Clayton copula \cr
		 \code{RotGumbel} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') 
		 }
}
\value{
CoVaR value evaluting at quantile CoVaR and quantile VaR.
}
\description{
Calculate the conditional quantile or CoVaR with different type of Copula and marginal distribution. In this package several bivariate copula families are included for bivariate analysis.
It provides functionality of elliptical (Gaussian and Student t) as well as Archimedean (Clayton, Gumbel, Frank, Plackett, BB1, SCJ, rotated clayton and rotated Gumbel) copulas to cover a large bandwidth of possible dependence structures.
}


\author{Andrea Ugolini <andreaugolini@me.com> \cr
        Juan Carlos Reboredo Noguiera <juancarlos.reboredo@usc.es>}

\references{Reboredo, J. C., & Ugolini, A. (2016). 
Quantile dependence of oil price movements and stock returns. Energy Economics, 54, 33-49.
}

\examples{
}

