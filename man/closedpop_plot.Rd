% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/closedpop_plot.R
\name{closedpop_plot}
\alias{closedpop_plot}
\title{A closed population ratio calculation function}
\usage{
closedpop_plot(maxage, Lmat, Lfish, M, Fi, Linf, k, a0, pW, qW, lambda,
  timeplot)
}
\arguments{
\item{maxage:}{max age of the species ie. number of age classes}

\item{Lmat:}{length at maturity}

\item{M:}{the natural mortality rate, if unknown generally use 0.2}

\item{Fi:}{the fishing mortality rate, F, find in stock assessment if don't have more localized estimate}

\item{Linf:}{asymptotic growth rate used in von-Bertallanfy growth equation}

\item{k:}{von-bertallanfy growth parameter estimate}

\item{a0:}{the age at length 0 used in the von-Bertallanfy growth equation}

\item{pW:}{weight length relationship parameter, same as a on fishbase.org but need to divide by 1000 to get in kg not grams}

\item{qW:}{weight length relationship parameter, same as b on fishbase.org}

\item{lambda:}{set the population growth rate in the MPA}

\item{timeplot:}{the time of the transient duration to plot on the output}
}
\description{
This function calculates the ratio change in a fished poulation after a marine protected area is implemented assuming  closed population
method is from White et al. 2013 'transient responses of fished populations to marine reserve establishment published in conservation letters
It uses a Leslies matrix and the output is a data frame with population ratio changes from a fished population to unfished in the MPA
}
\examples{
closedpop_plot(maxage=25,Lmat=18,Lfish=25,M=0.2,Fi=0.17, Linf=37.8,k=0.23,a0=-0.7,pW=6.29e-06,qW=3.172,lambda=1,timeplot=30)
}
\keyword{Leslie}
\keyword{closed}
\keyword{dynamics,}
\keyword{matrix}
\keyword{population}
\keyword{population,}
