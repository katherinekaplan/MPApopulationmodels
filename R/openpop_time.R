#' An open population function that calculates the length of time to reach equilibrium in response to an MPA
#'
#'This function calculates the ratio change in a fished poulation after a marine protected area is implemented assuming external recruitment.
#'The output is the time to reach the abundance and biomass equilibrium.
#' assumes a lambda of 1.0 in the MPA and calculate change from the stable age distribution of the fished population
#' @param maxage: max age of the species ie. number of age classes
#' @param Lmat: length at maturity
#' @param M: the natural mortality rate, if unknown generally use 0.2
#' @param Fi: the fishing mortality rate, F, find in stock assessment if don't have more localized estimate
#' @param Linf: von-bertallanfy growth parameter estimate, can find on fishbase
#' @param k: von-bertallanfy growth parameter estimate
#' @param a0: same as t0 in von-bertallanfy growth parameter
#' @param pw: weight length relationship estimate, same as a on fishbase.org but need to divide by 1000 to get in kg not grams
#' @param qw: weight length relationship estmate, same as b on fishbase.org
#' @param sig_r: stochastic parameter, log-normal distribution, around recruitment
#' @keywords closed population, population dynamics, Leslie matrix
#' @examples openpop_time(M=0.2,Fi=0.14,Lfish=25,Linf=37.8,k=0.13,a0=-0.7,maxage=25,pW=9.37e-06,qW=3.172)
#'openpop_time()


openpop_time = function(maxage,M,Fi,Lfish, Linf,k,a0,pW,qW) {
  tf=50
  R=500
  iterations=10
  MPAtime=5
  ##First step calculate the stable age distribution of the fished popultion
  a_harv0=(log((Lfish-Linf)/-Linf)/-k)+a0   ##age fished derived from length fished
  agefish=round(a_harv0,digits=0) ##make it a whole number to use in age structure matrix
  N0=rep(100,maxage) #Initial pop vector, start with 100 individual in each age class
  N0[1]=R ##except the first one, start with the atarting recruit number
  s=exp(-M)#no fishing case
  sf=exp(-(M+Fi)) ##fishing case
  sfx=rep(sf,maxage-1)
  sfx[1:agefish]=rep(s,agefish) ##set regular survival with no f for unfished age classes
  sxs=rep(s,maxage-1) #Survival vector ##number of s is ageclasses-1
  Nt = matrix(0,tf,maxage) #Initialize vector of population sizes with extra columns for spawners and recruitment before variability
  Nt[1,] = N0 #Put in initial values
  t<-1
  ##Get deterministic equilibrium
  for(t in 1:(tf-1)) {
    Nt[t+1,1] = R#*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
    Nt[t+1,2:(maxage)] = sfx*Nt[t,1:(maxage-1)] #Survivorship of each age class  in columns 2-10
  }
  ##Second step use that stable age disbribution of the fished population as the starting vector
  ##to determine the MPA effect
  N0=Nt[50,]##the stable age dist values from the fished state
  Nt2 = matrix(0,tf,maxage) #Initialize matrix of population sizes for second step
  Nt2[1,] = N0 #Put in initial values
  for(t in 1:(tf-1)) {
    ##For the first 5 time steps the population still fished
    if(t<=MPAtime){
      Nt2[t+1,1] = R#*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt2[t+1,2:(maxage)] = sfx*Nt2[t,1:(maxage-1)]
    }
    else{ ##then the population starts to fill in
      Nt2[t+1,1] = R#*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt2[t+1,2:(maxage)] = sxs*Nt2[t,1:(maxage-1)]#Survivorship of each age class  in columns 2-10
    }
  }
  final.N=rowSums(Nt2[,agefish:maxage]) ##include only fished age classes
  Nratio1=final.N/final.N[1]
  ##Now figure out the time point at which the final abundance is 95% of the equilibrium in the last time step
  ##Only works for no stochasticity
  final.N.ratio=final.N[tf]/final.N[1]
  time.ratio=final.N/final.N[tf]
  Ntime.to.equil=length(time.ratio[time.ratio<0.95])-5 #minus 5 for the first 5 years of running not in MPA
  Ntime.to.equil[Ntime.to.equil<0]=NA
  ##Next get biomass ratio change
  ##get lengths at age using the von-B equation
  a=seq(1:maxage)
  La=Linf*(1-exp(-k*(a-a0)))
  ##Now calculate weights at length
  w=pW*(La^qW)
  weights=Nt2[,agefish:maxage]%*%w[agefish:maxage]
  Bratio=weights/weights[1]
  final.B.ratio=weights[tf]/weights[1]
  time.ratio2=weights/weights[tf]
  Btime.to.equil=length(time.ratio2[time.ratio2<0.95])-5 #minus for the first 5 years of running not in MPA

  outputvars=data.frame(Ntime.to.equil,Btime.to.equil,final.N.ratio,final.B.ratio)
  colnames(outputvars)=c("Ntime.to.equil","Btime.to.equil","final.N.ratio","final.B.ratio")
  return( outputvars)
}



