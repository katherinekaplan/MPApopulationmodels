#' An open population ratio calculation function
#'
#'This function calculates the ratio change abundance and biomass of a fished poulation, after a marine protected area is implemented assuming a population with external recruitment.
#'The output is a data frame with population ratio changes from a fished population to unfished in the MPA
#'It includes output using stochastic recruitment
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
#' @examples openpop_ratio(M=0.2,Fi=0.14,Lfish=25,Linf=37.8,k=0.13,a0=-0.7,maxage=25,pW=9.37e-06,qW=3.172,sig_r=0.5)
#'openpop_ratio()

openpop_ratio = function(maxage,M,Fi,Lfish, Linf,k,a0,pW,qW,sig_r) {
  tf=50
   R=500
  iterations=10
  MPAtime=5
  ##First step calculate the stable age distribution of the fished popultion
  a_harv0=(log((Lfish-Linf)/-Linf)/-k)+a0   ##age fished
  agefish=round(a_harv0,digits=0)
  N0=rep(100,maxage) #Initial pop vector, start with 100 individual in each age class
  N0[1]=R
  #set.seed(1) #Set the seed so that every simulation uses same random sequence
  s=exp(-M)#no fishing case
  sf=exp(-(M+Fi)) ##fishing case
  sfx=rep(sf,maxage-1)
  sfx[1:agefish]=rep(s,agefish)
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
  ###Now add in the options for making several iterations with stochasticity
  ##make a function of within the function to work with and create iterations
  ##This starts the part for stochasticity, if iterations is 1 then no stochastic calcs
  ##First get a stable age dist with stochasticity to start the MPA response
  if (iterations>1){
    ##then create a function to do the iterations if that is specificed in the input
    AgeStructMatrix=function(tf,N0){
      ##get stabel age dist of fished state
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
        ##For the first 5 time steps the population still fished
        if(t<=5){
          Nt3[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
          Nt3[t+1,2:(maxage)] = sfx*Nt3[t,1:(maxage-1)] ##fishing case sfx for first 5 years
        }
        else{ ##then the population starts to fill in
          Nt3[t+1,1] = R*(exp(sig_r*rnorm(1,mean=0, sd=1))) #Recruits after variability in column 1, rnorm to generate random number for 1 point (n=1)
          Nt3[t+1,2:(maxage)] = sxs*Nt3[t,1:(maxage-1)]#Survivorship of each age class
        }}
      final.N3=rowSums(Nt3[,a_harv0:maxage]) ##include only fished age classes
      Nratio2=final.N3/final.N3[1]
      return(Nratio2)  ##returns your final response ratio
    }
    ## then repeat that function to get a distribution for the stochastic values
    itero=matrix(nrow=tf,ncol=iterations) ##this is the matrix to store the iterations
    i=1
    repeat{ ##this the repeat function
      itero[,i]=AgeStructMatrix(tf,N0)
      cat(itero[,i], "\n")
      i=i+1
      if(i>iterations) break
    }
    Nrat.mean=rowMeans(itero) ##this takes the mean of the iterations
    Nrat.sd=apply(itero,1,sd) ##this takes the standard deviation
    se=Nrat.sd/sqrt(iterations) ##this takes the standard error of the mean
    CI=2*se ##this is your confidence interval
    ###Now calculate stochasticity for weights
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
        ##For the first  time steps until MPAtime the population still fished
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
      weights=Nt4[,a_harv0:maxage]%*%w[a_harv0:maxage] ##get weights of fished age classes only
      bratio.stoch=weights/weights[1]
      return(bratio.stoch)
    }
    Brat.s=matrix(nrow=tf,ncol=iterations)
    i=1
    repeat{
      Brat.s[,i]=biomass.matrix(tf)
      cat(Brat.s[,i], "\n")
      i=i+1
      if(i>10) break
    }
    Brat.mean=rowMeans(Brat.s) ##this takes the mean for each iteration
    Brat.sd=apply(Brat.s,1,sd) ##this takes the standard deviation
    B.se=Brat.sd/sqrt(iterations) ##this takes the standard error of the mean
    B.CI=2*B.se
  }
  else { ##put NA in output if iterations is 1
    Nrat.mean=NA
    Nrat.sd=NA
    se=NA
    CI=NA
    Brat.mean=NA
    Brat.sd=NA
    B.se=NA
    B.CI=NA
  }


 Nratio=Nratio1
  Nrat.itero.mean=Nrat.mean
  Nrat.itero.CI=CI
   Bratio.itero.mean=Brat.mean
 Bratio.itero.CI= B.CI
  time=seq(1,tf)
  outputvars=data.frame(time,Nratio,Bratio,Nrat.itero.mean,Nrat.itero.CI,Bratio.itero.mean,
                        Bratio.itero.CI)
  colnames(outputvars)=c("time","Nratio","Bratio","Nrat.itero.mean","Nrat.itero.CI","Bratio.itero.mean","Bratio.itero.CI")

  return(outputvars)
}




