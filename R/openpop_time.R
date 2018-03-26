#' An open population function that calculates the length of time to reach abundance and biomass equilibrium in response to a marine protected area, MPA
#'
#'This function calculates the ratio change in a fished population after a marine protected area is implemented assuming external recruitment.
#'The output is the approximate time to reach the abundance and biomass equilibrium. For long-lived species used longer timeframes to ensure final ratios are reached.
#' @param tf time steps to run the population
#' @param maxage: max age of the species ie. number of age classes
#' @param Lmat: length at maturity
#' @param M: the natural mortality rate, if unknown generally use 0.2
#' @param Fi: the fishing mortality rate, F, find in stock assessment if don't have more localized estimate
#' @param Linf: asymptotic length of the fish
#' @param k: von-bertallanfy growth parameter estimate
#' @param a0: the age at length 0 used in the von-Bertallanfy growth equation
#' @param pW: weight length relationship parameter, same as a on fishbase.org but need to divide by 1000 to get in kg not grams
#' @param qW: weight length relationship parameter, same as b on fishbase.org
#' @param R: the number of recruits entering the population
#' @return N.time.to.equil= the number of years it will take for the population abundance to reach its final equilibrium abundance
#' @return Btime.to.equil= the number of years it will take for the population biomass to reach its final equilibirum
#' @return final.N.ratio= the final abundance ratio increase reached in comparison to pre-MPA abundance for fished age classes
#' @return final.B.ratio= the final biomass ratio increase reach in compared to pre-MPA biomass for fished age classes
#' @keywords open population, population dynamics
#' @examples openpop_time(tf=50,M=0.2,Fi=0.14,Lfish=25,Linf=37.8,k=0.13,a0=-0.7,maxage=25,pW=9.37e-06,qW=3.172,R=500)

openpop_time = function(tf,maxage,M,Fi,Lfish, Linf,k,a0,pW,qW,R) {
  ##First step calculate the stable age distribution of the fished popultion
  a_harv0=(log((Lfish-Linf)/-Linf)/-k)+a0   ##age fished back calculated from von-B eqn
  agefish=round(a_harv0,digits=0) ##round age fished to integer
  N0=rep(100,maxage) #Initial pop vector, start with 100 individual in each age class
  N0[1]=R
  s=exp(-M)#no fishing case
  sf=exp(-(M+Fi)) ##fishing case
  sfx=rep(sf,maxage-1)
  sfx[1:(agefish-1)]=rep(s,(agefish-1))
  sxs=rep(s,maxage-1) #Survival vector ##number of s is ageclasses-1
  Nt = matrix(0,tf,maxage) #Initialize vector of population sizes with extra columns for spawners and recruitment before variability
  Nt[1,] = N0 #Put in initial values
  set.seed(1) #Set the seed so that every simulation uses same random sequence
  t<-1
  ##Get deterministic equilibrium for fished state
  for(t in 1:(tf-1)) {
    Nt[t+1,1] = R#*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
    Nt[t+1,2:(maxage)] = sfx*Nt[t,1:(maxage-1)] #Survivorship of each age class  in columns 2-10
  }
  ##Second step use that stable age disbribution of the fished population as the starting vector
  ##to determine the MPA effect
  N0=Nt[tf,]##the stable age dist values from the fished state
  N.fished=N0
  Nt2 = matrix(0,tf,maxage) #Initialize matrix of population sizes for second step
  Nt2[1,] = N0 #Put in initial values
  for(t in 1:(tf-1)) {
    MPAtime=0
    if(t<=MPAtime){
      Nt2[t+1,1] = R#*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt2[t+1,2:(maxage)] = sfx*Nt2[t,1:(maxage-1)]
    }
    else{ ##then the population starts to fill in
      Nt2[t+1,1] = R#*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt2[t+1,2:(maxage)] = sxs*Nt2[t,1:(maxage-1)]#Survivorship of each age class  in columns 2-10
    }
  }
  Nt2=Nt2
  final.N=rowSums(Nt2[,agefish:maxage]) ##include only fished age classes
  Nratio1=final.N/final.N[1]
  ##Now figure out the time point at which the final abundance is 95% of the equilibrium in the last time step
  ##Only works for no stochasticity
  final.N.ratio=final.N[tf]/final.N[1]
  time.ratio=final.N/final.N[tf]
  Ntime.to.equil=length(time.ratio[time.ratio<0.95])-MPAtime #minus the time when MPA implemented
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
  Btime.to.equil=length(time.ratio2[time.ratio2<0.95])-MPAtime #minus time to MPa implementation of running not in MPA

  outputvars=data.frame(Ntime.to.equil,Btime.to.equil,final.N.ratio,final.B.ratio)
  colnames(outputvars)=c("Ntime.to.equil","Btime.to.equil","final.N.ratio","final.B.ratio")
  return( outputvars)
}


