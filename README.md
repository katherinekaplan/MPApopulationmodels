# MPApopulationmodels

The R package 'MPApopulationmodels' contains several functions that can be used to project population responses to marine protected areas (MPAs).  Input parameters including life history characteristics and a prior fishing mortality rate are required to use the functions described. Each function is also documented in R and can be searched once the package is downloaded into R by typing ??function_name. The input and output for each function can be found in the R documentation.  Below I provide a brief summary of the utility of each function. 

The following functions are available from this package:

closedpop_metrics: calculates the ratio change in a fished population after a marine protected area is implented assuming a closed population. It provides P1, the period of oscillations; rho, the rate of return to the stable age distribution in the MPA; theta, the angle between the fished state population vector and the stable age distribution in the MPA; Nratio, the abundance changes over time; Bratio, the biomass ratio change over time; and transient_length, the approximate length of the transient period calculated from rho.

closedpop_plot: creates a plot using ggplot that shows the changes in abundance and biomass expected following MPA implementation

closedpop_ratio: returns a data frame with time, abundance ratio (Nratio), and biomass ratio (Bratio), which is used by the closedpop_plot function to generate the plot shown

openpop_plots: returns 3 plots, the first is the deterministic population projection of abundance and biomass changing over time, the second plot is a population projection with stochastic recruitment and returns the upper and lower quartile for changes in abundance overtime since MPA implementation based on the number of simulations specified in the input of the function, the third plot is a stochastic population projection also but calculates changes in biomass over time.

openpop_ratio: outputs a data frame that contains changes in the population abundance and biomass overtime since MPA implementation as shown with the openpop_plots function.  The data frame includes a column for time, deterministic abundance ratio changes (Nratio), deterministic biomass ratio changes (Bratio), the mean of simulations specified for abundance changes with stochasticity in recruitment inside an MPA (Nrat.sim.mean), the lower quartile for the distribution of simulated runs for changes in abundance in an MPA with stochasticity in recruitment (Nrat.lowerCI.MPA), the upper quartile for the distribution of simulated runs for changes in abundance in an MPA with stochasticity in recruitment (Nrat.upperCI.MPA), the mean of a simulations of abundance changes with stochasticity in recruitment for a control (fished site) not inside an MPA (Nrat.mean.noMPA), the lower quartile for the distribution of simulated runs for changes in abundance in a control (fished site) with stochasticity in recruitment (Nrat.lowerCI.noMPA), the upper quartile for the distribution of simulated runs for changes in abundance in a control (fished site) with stochasticity in recruitment (Nrat.upperCI.noMPA),the mean of a simulations of biomass changes with stochasticity in recruitment inside an MPA (Bratio.sim.mean), the lower quartile for the distribution of simulated runs for changes in biomass in an MPA with stochasticity in recruitment (Brat.lowerCI.MPA), the upper quartile for the distribution of simulated runs for changes in biomass in an MPA with stochasticity in recruitment (Brat.upperCI.MPA), the mean of simulations specified for biomass changes with stochasticity in recruitment inside an MPA (Bratio.sim.mean),the mean of simulations specified for biomass changes with stochasticity in recruitment in a control site (Brat.mean.noMPA),the lower quartile for the distribution of simulated runs for changes in biomass in a control site with stochasticity in recruitment (Brat.lowerCI.noMPA), and the upper quartile for the distribution of simulated runs for changes in biomass in an MPA with stochasticity in recruitment (Brat.upperCI.noMPA).

openpop_time: provides outputs based on determinstic population projections based on the time reach 95% of the final equilibrium ratio in response to an MPA.  Outputs include the abundance time to equilbirium (Ntime.to.equil), the biomass time to equilibrium (Btime.to.equil), the final abundance ratio change (final.N.ratio), and the final biomass ratio change (final.B.ratio).
