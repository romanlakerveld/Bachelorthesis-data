R Notebook
This is the first model by Martinez et al. displaying the interactions

library("ggplot2")
library("dplyr")
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library("deSolve")
library('tidyr')
library('directlabels')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Martinez1 <- function(t, y, parameters) {
  #calculates dx/dts for a GRN
  
  # t: time at which to evaluate derivatives
  # y: vector of system variables (c(P,B,R))
  # parameters: vector of model parameters
  
  BLIMP1 <- y[1]
  BCL6 <- y[2]
  IRF4 <- y[3]
  
  up <- parameters['up']    # passive transcription rate
  ub <- parameters['ub']
  ur <- parameters['ur']
  
  op <- parameters['op']    # max induced transcription rate
  ob <- parameters['ob']
  or <- parameters['or']
  
  kb <- parameters['kb']    # dissociation constant
  kr <- parameters['kr']
  kp <- parameters['kp']
  
  ep <- parameters['ep']    # rate of degradation
  eb <- parameters['eb']
  er <- parameters['er']
  
  CD40 <- dnorm(t, 30, 2) * (kb ^ 2 / (kb ^ 2 + BCL6 ^ 2))
  BCR <- dnorm(t, 40, 5) * (kb ^ 2 / ( kb ^ 2 + BCL6 ^ 2))
  
  # calculate rate of change
  dBLIMP1 <- up + op * (kb ^ 2 / (kb ^ 2 + BCL6 ^ 2)) + op * (IRF4 ^ 2 / (kr ^ 2 + IRF4 ^ 2)) - ep * BLIMP1
  dBCL6 <- ub + ob * (kp ^ 2 / (kp ^ 2 + BLIMP1 ^ 2)) * (kb ^ 2 / (kb ^ 2 + BCL6 ^ 2)) * (kr ^ 2 / (kr ^ 2 + IRF4 ^ 2)) - (eb + BCR) * BCL6
  dIRF4 <- ur + or * (IRF4 ^ 2 / (kr ^ 2 + IRF4 ^ 2)) + CD40 - er * IRF4
  
  # return rate of change
  return(list(c(dBLIMP1, dBCL6, dIRF4)))
}

# run the numerical solution

parameters = c(
  up = 10 ^ -6,
  ub = 2,
  ur = 0.1,
  op = 9,
  ob = 100,       #just some parameters, no bxiggie
  or = 2.6,
  kp = 1,
  kb = 1,
  kr = 1,
  ep = 1,
  eb = 1,
  er = 1,
  BCR = 0,
  CD40 = 0
)

state <- c(BLIMP1 = 0.5, BCL6 = 5, IRF4 = 0.2) # starting states

times <- seq(0, 100, by = 0.01)

result <- ode(y=state, times = times, func = Martinez1, parms = parameters)
result <- data.frame(result)

result <- mutate(result, BCR = dnorm(time, 40, 5) * 150 * (1 / ( 1 + BCL6 ^ 2))) %>%
  mutate(CD40 = dnorm(time, 50, 5) * (1 / ( 1 + BCL6 ^ 2)))

result <- result %>% 
  gather(Column, Value, -time) %>% 
  filter(Column != "BCR" & Column != "CD40")

# plot the results

result %>% 
  ggplot(aes(x=time, y= Value, color = Column)) +
  scale_colour_manual(values=cbPalette) +
  geom_line(aes(x = time, y = Value)) +
  scale_fill_hue(l=20) +
  labs(y = "Concentration") +
  geom_dl(aes(label = Column), method = list(dl.combine("last.points"), cex = 0.8))
