#' An open population function to calculate population changes in a marine protected area, MPA
#'
#'This function calculates the ratio change abundance and biomass of a fished poulation, after a marine protected area is implemented assuming a population with external recruitment.
#'The output is a data frame with population ratio changes from a fished population to unfished in the MPA
#'It includes deterministic population model output and output using stochastic recruitment
#' @param tf: the time steps to run the population
#' @param maxage: max age of the species ie. number of age classes
#' @param Lmat: length at maturity
#' @param M: the natural mortality rate
#' @param Fi: the fishing mortality rate, F, find in stock assessment if don't have more localized estimate
#' @param Linf: von-bertallanfy growth parameter estimate, can find on fishbase
#' @param k: von-bertallanfy growth parameter estimate
#' @param a0: same as t0 in von-bertallanfy growth parameter
#' @param pw: weight length relationship estimate, same as a on fishbase.org but need to divide by 1000 to get in kg not grams
#' @param qw: weight length relationship estmate, same as b on fishbase.org
#' @param R: number of recruits entering the population
#' @param sig_r: stochastic parameter, log-normal distribution, around recruitment
#' @param MPAtime: the time step to implement the MPA
#' @param simulations: the number of simulations to run
#' @return Nratio: the abundance ratio change over time
#' @return Bratio: the biomass ratio change over time
#' @return Nrat.sim.mean: the mean of abundance for simulation runs with stochastic recruitment in the MPA
#' @return Nrat.lowerCI.MPA: the lower quartile of abundance runs with stochastic recruitment in the MPA
#' @return Nrat.upperCI.MPA: the upper quartile of abundance runs with stochastic recruitment in the MPA
#' @return Nrat.mean.noMPA: the mean of changes in simulated abundance with stochastic recruitment in the fished state
#' @return Nrat.lowerCI.noMPA: the lower quartile of simulated runs for abundance with stochastic recruitment in the fished state
#' @return Nrat.upperCI.noMPA: the lower quartile of simulated runs with stochastic recruitment in the fished state
#' @return Bratio.sim.mean: the mean of changes in biomass for simulation runs with stochastic recruitment in the MPA
#' @return Brat.lowerCI.MPA: the lower quartile of biomass runs with stochastic recruitment in the MPA
#' @return Brat.upperCI.MPA: the upper quartile of biomass runs with stochastic recruitment in the MPA
#' @return Brat.mean.noMPA: the mean of changes in biomass for simulation runs with stochastic recruitment in the fished state
#' @return Brat.lowerCI.noMPA: the lower quartile of biomass runs with stochastic recruitment in the fished state
#' @return Brat.upperCI.noMPA: the upper quartile of biomass runs with stochastic recruitment in the fished state
#' @keywords open population, population dynamics
#' @examples openpop_ratio(tf=50, M=0.2,Fi=0.14,Lfish=25,Linf=37.8,k=0.13,a0=-0.7,maxage=25,pW=9.37e-06,qW=3.172,R=500,
#'  sig_r=0.5, MPAtime=1,simulations=200)


openpop_ratio = function(tf,maxage,M,Fi,Lfish, Linf,k,a0,pW,qW,R,sig_r,MPAtime,simulations) {
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
  ###Now add in the options for making several simulations with stochasticity
  ##This starts the part for stochasticity, if simulations is >1 then no stochastic calcs
  ##First get a stable age dist with stochasticity to start the MPA response
  ##then create a function to do the simulations if that is specificed in the input
  AgeStructMatrix=function(tf,N0){
    ##get stable age dist of fished state
    Nt0 = matrix(0,tf,maxage) #Initialize vector of population sizes
    Nt0[1,] = rep(100,maxage) #Put in initial values just a starting vector
    Nt0[1,1]=500
    for(t in 1:(tf-1)) {
      Nt0[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability  rnorm to generate random number for 1 point (n=1)
      Nt0[t+1,2:(maxage)] = sfx*Nt0[t,1:(maxage-1)] #Survivorship of each age class
    }
    N0=Nt0[tf,]##the stable age dist values from the fished state again, but now with stochasticity
    Nt3 = matrix(0,tf,maxage) #Initialize matrix of population sizes for second step
    Nt3[1,] = N0 #Put in initial values
    for(t in 1:(tf-1)) {
      ##For the first time steps until MPA specified the population still fished
      if(t<=MPAtime){
        Nt3[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
        Nt3[t+1,2:(maxage)] = sfx*Nt3[t,1:(maxage-1)] ##fishing case sfx for first 5 years
      }
      else{ ##then the population starts to fill in
        Nt3[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
        Nt3[t+1,2:(maxage)] = sxs*Nt3[t,1:(maxage-1)]#Survivorship of each age class
      }}
    final.N3=rowSums(Nt3[,a_harv0:maxage]) ##include only fished age classes
    Nratio2=final.N3/sum(Nt0[tf-1,a_harv0:maxage])
    return(Nratio2)  ##returns your final response ratio
  }
  ## then repeat that function to get a distribution for the stochastic values
  itero=matrix(nrow=tf,ncol=simulations) ##this is the matrix to store the simulations
  i=1
  repeat{ ##this the repeat function
    itero[,i]=AgeStructMatrix(tf,N0)
    cat(itero[,i], "\n")
    i=i+1
    if(i>simulations) break
  }
  Nrat.mean=rowMeans(itero) ##this takes the mean of the simulations
  Nrat.quantiles.MPA=apply(itero,1,quantile,probs=c(0.25,0.75))
  NoMPAMatrix=function(tf,N0){
    ##if there were no MPA case for plot
    ##first run with random recruitment for 50 years
    Nt0 = matrix(0,tf,maxage) #Initialize vector of population sizes
    for(t in 1:(tf-1)) {
      Nt0[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt0[t+1,2:(maxage)] = sfx*Nt0[t,1:(maxage-1)] #Survivorship of each age class  in columns 2-10
    }
    N0=Nt0[tf,]##the stable age dist values from the fished state again, but now with stochasticity
    Nt.noMPA = matrix(0,tf,maxage) #Initialize matrix of population sizes for second step
    Nt.noMPA[1,] = N0 #Put in initial values
    for(t in 1:(tf-1)) {
      Nt.noMPA[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability  rnorm to generate random number for 1 point (n=1)
      Nt.noMPA[t+1,2:(maxage)] = sfx*Nt.noMPA[t,1:(maxage-1)] #Survivorship of each age class
    }
    final.Nt.noMPA=rowSums(Nt.noMPA[,a_harv0:maxage]) ##include only fished age classes
    Nratio.noMPA=final.Nt.noMPA/sum(Nt0[tf-1,a_harv0:maxage])
    return(Nratio.noMPA)
  }

  ## then repeat that function to get a distribution for the stochastic values
  itero2=matrix(nrow=tf,ncol=simulations) ##this is the matrix to store the simulations
  i=1
  repeat{ ##this the repeat function
    itero2[,i]= NoMPAMatrix(tf,N0)
    cat(itero2[,i], "\n")
    i=i+1
    if(i>simulations) break
  }
  Nrat.mean.noMPA=rowMeans(itero2) ##this takes the mean of the simulations
  Nrat.quantiles.noMPA=apply(itero2,1,quantile,probs=c(0.25,0.75))
  ###Do the same thing for biomass
  ###Now calculate stochasticity for weights with MPA
  ##first get weights
  biomass.matrix=function(tf){
    ##Now get starting abundance of fished state
    Nt0 = matrix(0,tf,maxage) #Initialize vector of population sizes
    Nt0[1,] = rep(100,maxage) #Put in initial values just a starting vector
    for(t in 1:(tf-1)) {
      Nt0[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt0[t+1,2:(maxage)] = sfx*Nt0[t,1:(maxage-1)] #Survivorship of each age class
    }
    Nt4 = matrix(0,tf,maxage)
    Nt4[1,] = Nt0[tf,]##starting abundance vector from the fished state
    for(t in 1:(tf-1)) {
      ##For the first time steps the population still fished until MPA starts
      if(t<=MPAtime){
        Nt4[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
        Nt4[t+1,2:(maxage)] = sfx*Nt4[t,1:(maxage-1)]
      }
      else{ ##then the population starts to fill in
        Nt4[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
        Nt4[t+1,2:(maxage)] = sxs*Nt4[t,1:(maxage-1)]#Survivorship of each age class  in columns 2-10
      }}
    a=seq(1:maxage)
    La=Linf*(1-exp(-k*(a-a0))) ##von-B eqn calculates lengths at age
    ##Now calculate weights at length
    w=pW*(La^qW) ##weights at length in kg
    weight1=Nt0[tf-1,a_harv0:maxage]%*%w[a_harv0:maxage]
    weights=Nt4[,a_harv0:maxage]%*%w[a_harv0:maxage] ##get weights of fished age classes only
    bratio.stoch=weights/weight1[1]
    return(bratio.stoch)
  }
  itero3=matrix(nrow=tf,ncol=simulations)
  i=1
  repeat{
    itero3[,i]=biomass.matrix(tf)
    cat(itero3[,i], "\n")
    i=i+1
    if(i>simulations) break
  }
  Brat.mean=rowMeans(itero3) ##this takes the mean for each iteration
  Brat.quantiles.MPA=apply(itero3,1,quantile,probs=c(0.25,0.75))

  ##Now biomass with no MPA
  NoMPAMatrix.biomass=function(tf,N0){
    ##if there were no MPA case for plot
    ##first run with random recruitment for 50 years
    Nt0 = matrix(0,tf,maxage) #Initialize vector of population sizes
    for(t in 1:(tf-1)) {
      Nt0[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
      Nt0[t+1,2:(maxage)] = sfx*Nt0[t,1:(maxage-1)] #Survivorship of each age class
    }
    N0=Nt0[tf,]##the stable age dist values from the fished state again, but now with stochasticity
    Nt.noMPA = matrix(0,tf,maxage) #Initialize matrix of population sizes for second step
    Nt.noMPA[1,] = N0 #Put in initial values
    for(t in 1:(tf-1)) {
      Nt.noMPA[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability  rnorm to generate random number for 1 point (n=1)
      Nt.noMPA[t+1,2:(maxage)] = sfx*Nt.noMPA[t,1:(maxage-1)] #Survivorship of each age class
    }
    a=seq(1:maxage)
    La=Linf*(1-exp(-k*(a-a0))) ##von-B eqn calculates lengths at age
    ##Now calculate weights at length
    w=pW*(La^qW) ##weights at length in kg
    weights2=Nt.noMPA[,a_harv0:maxage]%*%w[a_harv0:maxage] ##get weights of fished age classes only
    weight1=Nt0[tf-1,a_harv0:maxage]%*%w[a_harv0:maxage]
    bratio.stoch.noMPA=weights2/weight1[1]
    return(bratio.stoch.noMPA)
  }
  ## then repeat that function to get a distribution for the stochastic values
  itero4=matrix(nrow=tf,ncol=simulations) ##this is the matrix to store the simulations
  i=1
  repeat{ ##this the repeat function
    itero4[,i]= NoMPAMatrix.biomass(tf,N0)
    cat(itero4[,i], "\n")
    i=i+1
    if(i>simulations) break
  }
  Brat.mean.noMPA=rowMeans(itero4) ##this takes the mean of the simulations
  Brat.quantiles.noMPA=apply(itero4,1,quantile,probs=c(0.25,0.75))
  ##put outputs in data frame to plot
  t=seq(1,tf)
  newdf=data.frame(t,Nratio1,Bratio,
                   Nrat.mean,Nrat.quantiles.MPA[1,],Nrat.quantiles.MPA[2,],
                   Nrat.mean.noMPA,Nrat.quantiles.noMPA[1,],Nrat.quantiles.noMPA[2,],
                   Brat.mean,Brat.quantiles.MPA[1,], Brat.quantiles.MPA[2,],
                   Brat.mean.noMPA, Brat.quantiles.noMPA[1,],Brat.quantiles.noMPA[2,])
  colnames(newdf)=c("time","Nratio","Bratio",
                    "Nrat.sim.mean","Nrat.lowerCI.MPA","Nrat.upperCI.MPA",
                    "Nrat.mean.noMPA","Nrat.lowerCI.noMPA","Nrat.upperCI.noMPA",
                    "Bratio.sim.mean","Brat.lowerCI.MPA","Brat.upperCI.MPA",
                    "Brat.mean.noMPA","Brat.lowerCI.noMPA","Brat.upperCI.noMPA")
  return(newdf)
}




