#' A closed population ratio calculation function
#'
#'This function calculates the ratio change in a fished poulation after a marine protected area is implemented assuming  closed population
#'method is from White et al. 2013 'transient responses of fished populations to marine reserve establishment published in conservation letters
#'It uses a Leslies matrix and the output is a data frame with population ratio changes from a fished population to unfished in the MPA

#' @param maxage: max age of the species ie. number of age classes
#' @param Lmat: length at maturity
#' @param M: the natural mortality rate, if unknown generally use 0.2
#' @param Fi: the fishing mortality rate, F, find in stock assessment if don't have more localized estimate
#' @param Linf: von-bertallanfy growth parameter estimate, can find on fishbase
#' @param k: von-bertallanfy growth parameter estimate
#' @param a0: same as t0 in von-bertallanfy growth parameter
#' @param pw: weight length relationship estimate, same as a on fishbase.org but need to divide by 1000 to get in kg not grams
#' @param qw: weight length relationship estmate, same as b on fishbase.org
#' @keywords closed population, population dynamics, Leslie matrix
#' @examples closedpop_ratio(maxage=25,Lmat=18,Lfish=25,M=0.2,Fi=0.17, Linf=37.8,k=0.23,a0=-0.7,pW=6.29e-06,qW=3.172)

closedpop_ratio = function(maxage,Lmat,Lfish,M,Fi, Linf,k,a0,pW,qW) {
  tf=50 #time steps to run the population
  time2=seq(1,50)
  a_mat0=(log((Lmat-Linf)/-Linf)/-k)+a0 ##calculate the age at maturity from length
  a_mat=round(a_mat0,digits=0)##rounds that age to whole number to input into leslie matrix
  a_harv0=(log((Lfish-Linf)/-Linf)/-k)+a0   ##age fished calculated from length fished
  a_harv=round(a_harv0,digits=0)##rounds to whole number to input into leslie matrix
  amf=a_mat0/a_harv0 ##ratio of age at maturity to age fished
  fm_ratio=Fi/M ##ratio fo fishing mortality to natural mortality
  # repro and mortality
  Ac = c(rep(0,a_mat),rep(1,(maxage-a_mat)))       # set age at first reproduction vector
  M = c(rep(M,(maxage)),0)		               # vector of natural mortality at age a
  H = c(rep(0,a_harv),rep(Fi,maxage-a_harv)) # vector of fishing mortality at age a
  ma = array(0,dim=c(maxage))                    # set base array for fecundity at age
  La = ma                                    # " length at age
  Wa = ma                                    # " weight at age
  for (a in 1:(maxage)) {
    La[a] = Linf*(1-exp(-k*(a-a0)))                # vonBert growth function
    Wa[a] = pW*(La[a])^qW                           # weight-length relationship  # Weight-length coefficient kg*cm^q
    ma[a] = Wa[a]*Ac[a]                             # fecundity proportional to weight
  }
  ##Alpha calc function, calculate the alpha (fecundity value)
  ##needed to keep the unfished population growing at at lambda rate of 1.02
  ##can adjust that growth rate assumption by changing lambdaGoal
  lambdaGoal=1.0
  alphafun=function(alpha){ ##function to calculate lambda, which will be used below to calculate correct alpha level
    ma_final = alpha*ma ##top row of Leslie matrix
    projM_unfished = rbind(ma_final,cbind(diag(c(exp(-M[1:(maxage-1)]))),0)) ##sets up Leslie matrix
    lambda = Re(eigen(projM_unfished)$values[1]) ##Lambda calculation is 1st dominant eigenvalue
    return(lambda)
  }

  alphaMin=function(alpha) {  ## minimize the alphafun above to find appropriate alpha level
    err=(lambdaGoal-alphafun(alpha)) ##error is the difference between the goal lambda and our function of lambda
    sum(err^2)/length(err) ##sum of squared error to minimize in optimization function
  }

  OptimOut=optimize(alphaMin,lower=0,upper=5) ##uses optimization function to determine alpha level
  alpha=OptimOut$minimum
  ##gives correct alpha value
  ma_final = alpha*ma                            # final fecundity values to go across top of Leslie matrix
  # construct Leslie matrix by combining ma with diagonal matrix of survival
  projM_unfished = rbind(ma_final,cbind(diag(c(exp(-M[1:(maxage-1)]))),0))  # diag creates diagonal matrix with terms exp(-M) for a = 1:maxage-1 and 0 for a = maxage
  projM_fished = rbind(ma_final,cbind(diag(c(exp(-(M[1:(maxage-1)]+H[1:(maxage-1)])))),0)) #MPA Leslie matrix
  lambda1.unfished = eigen(projM_unfished)$values[1]                 # lambda1 MPA
  lambda2.unfished = eigen(projM_unfished)$values[2]                 # lambda2 MPA
  lambda1.fished = eigen(projM_fished)$values[1]                # lambda1 fished
  lambda2.fished = eigen(projM_fished)$values[2]                 # lambda2 fished
  # calculate 1st right eigenvector of projM to get SAD, v1
  ##v1 is the fished stable age distribution
  v1 = abs((Re(eigen(projM_fished)$vectors[,1])))    # real portion of eigenvector
  N0=rep(1000,maxage)
  N0=N0*v1 ##get starting pop vector by multiplying by the fished stable age dist
  Nc = matrix(0,tf,maxage) #Initialize vector of population sizes
  Nc[1,] = N0 #Put in initial values
  for(t in 1:(tf-1)) {
    Nc[t+1,1:(maxage)] =projM_unfished%*%Nc[t,]#calculates population changing over time with MPA
  }
  final.Nc=rowSums(Nc[,a_harv:maxage]) ##include only fished age classes in final abundance ratios
  Nratioc=final.Nc/final.Nc[1]
  ##Calculate weights at length
  ##Now calculate weights at length
  w=pW*(La^qW)
  weights=Nc[,a_harv:maxage]%*%w[a_harv:maxage]
  Bratio=weights/weights[1]
  ##Calculate period
  P1=2*pi/atan((Im(lambda2.fished)/Re(lambda2.fished)))
  ##Calculate the rate of convergence, rho
  rho=Re(lambda1.fished/abs(lambda2.fished))
  ##calculate total length of transience
  t.t=log(0.1)/-log(rho)
  ##Calculate theta
  w1=abs(Re(eigen(projM_unfished)$vectors[,1]))
  theta=acos((N0%*%w1)/(sqrt(sum(N0^2))*sqrt(sum(w1^2))))
  ages=seq(1:maxage)
  mean.spawn=sum(((v1*ma_final)/sum(v1*ma_final))*ages)
  Cv.spawn=sd(((v1*ma_final)/sum(v1*ma_final))*ages)/mean.spawn
  output=data.frame(Nratioc=Nratioc[1:30],time2=time2[1:30])
  return(output)
}


