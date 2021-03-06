% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/openpop_time.R
\name{openpop_time}
\alias{openpop_time}
\title{An open population function that calculates the length of time to reach abundance and biomass equilibrium in response to a marine protected area, MPA}
\usage{
openpop_time(tf, maxage, M, Fi, Lfish, Linf, k, a0, pW, qW, R)
}
\arguments{
\item{tf}{time steps to run the population}

\item{maxage:}{max age of the species ie. number of age classes}

\item{Lmat:}{length at maturity}

\item{M:}{the natural mortality rate, if unknown generally use 0.2}

\item{Fi:}{the fishing mortality rate, F, find in stock assessment if don't have more localized estimate}

\item{Linf:}{asymptotic length of the fish}

\item{k:}{von-bertallanfy growth parameter estimate}

\item{a0:}{the age at length 0 used in the von-Bertallanfy growth equation}

\item{pW:}{weight length relationship parameter, same as a on fishbase.org but need to divide by 1000 to get in kg not grams}

\item{qW:}{weight length relationship parameter, same as b on fishbase.org}

\item{R:}{the number of recruits entering the population}
}
\value{
N.time.to.equil= the number of years it will take for the population abundance to reach its final equilibrium abundance

Btime.to.equil= the number of years it will take for the population biomass to reach its final equilibirum

final.N.ratio= the final abundance ratio increase reached in comparison to pre-MPA abundance for fished age classes

final.B.ratio= the final biomass ratio increase reach in compared to pre-MPA biomass for fished age classes
}
\description{
This function calculates the ratio change in a fished population after a marine protected area is implemented assuming external recruitment.
The output is the approximate time to reach the abundance and biomass equilibrium. For long-lived species used longer timeframes to ensure final ratios are reached.
}
\examples{
openpop_time(tf=50,M=0.2,Fi=0.14,Lfish=25,Linf=37.8,k=0.13,a0=-0.7,maxage=25,pW=9.37e-06,qW=3.172,R=500)
}
\keyword{dynamics}
\keyword{open}
\keyword{population}
\keyword{population,}
