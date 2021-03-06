% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/openpop_ratio.R
\name{openpop_ratio}
\alias{openpop_ratio}
\title{An open population function to calculate population changes in a marine protected area, MPA}
\usage{
openpop_ratio(tf, maxage, M, Fi, Lfish, Linf, k, a0, pW, qW, R, sig_r, MPAtime,
  simulations)
}
\arguments{
\item{tf:}{the time steps to run the population}

\item{maxage:}{max age of the species ie. number of age classes}

\item{Lmat:}{length at maturity}

\item{M:}{the natural mortality rate}

\item{Fi:}{the fishing mortality rate, F, find in stock assessment if don't have more localized estimate}

\item{Linf:}{asymptotic growth rate used in von-Bertallanfy growth equation, can find on fishbase}

\item{k:}{von-bertallanfy growth parameter estimate}

\item{a0:}{the age at length 0 used in the von-Bertallanfy growth equation}

\item{pW:}{weight length relationship parameter, same as a on fishbase.org but need to divide by 1000 to get in kg not grams}

\item{qW:}{weight length relationship parameter, same as b on fishbase.org}

\item{R:}{number of recruits entering the population}

\item{sig_r:}{stochastic parameter, log-normal distribution, around recruitment}

\item{MPAtime:}{the time step to implement the MPA}

\item{simulations:}{the number of simulations to run}
}
\value{
Nratio: the abundance ratio change over time

Bratio: the biomass ratio change over time

Nrat.sim.mean: the mean of abundance for simulation runs with stochastic recruitment in the MPA

Nrat.lowerCI.MPA: the lower quartile of abundance runs with stochastic recruitment in the MPA

Nrat.upperCI.MPA: the upper quartile of abundance runs with stochastic recruitment in the MPA

Nrat.mean.noMPA: the mean of changes in simulated abundance with stochastic recruitment in the fished state

Nrat.lowerCI.noMPA: the lower quartile of simulated runs for abundance with stochastic recruitment in the fished state

Nrat.upperCI.noMPA: the lower quartile of simulated runs with stochastic recruitment in the fished state

Bratio.sim.mean: the mean of changes in biomass for simulation runs with stochastic recruitment in the MPA

Brat.lowerCI.MPA: the lower quartile of biomass runs with stochastic recruitment in the MPA

Brat.upperCI.MPA: the upper quartile of biomass runs with stochastic recruitment in the MPA

Brat.mean.noMPA: the mean of changes in biomass for simulation runs with stochastic recruitment in the fished state

Brat.lowerCI.noMPA: the lower quartile of biomass runs with stochastic recruitment in the fished state

Brat.upperCI.noMPA: the upper quartile of biomass runs with stochastic recruitment in the fished state
}
\description{
This function calculates the ratio change abundance and biomass of a fished poulation, after a marine protected area is implemented assuming a population with external recruitment.
The output is a data frame with population ratio changes from a fished population to unfished in the MPA
It includes deterministic population model output and output using stochastic recruitment
}
\examples{
openpop_ratio(tf=50, M=0.2,Fi=0.14,Lfish=25,Linf=37.8,k=0.13,a0=-0.7,maxage=25,pW=9.37e-06,qW=3.172,R=500,
 sig_r=0.5, MPAtime=1,simulations=100)
}
\keyword{dynamics}
\keyword{open}
\keyword{population}
\keyword{population,}
