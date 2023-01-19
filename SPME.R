# Load libraries
library(dplyr) # organize data
library(tidyverse)
library(ggplot2) # plotting
library(deSolve) # solving differential equations

# Reactive transport function ---------------------------------------------

rtm.SPME = function(t, c, parms){
  
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  
  # Constant passed through a list called parms
  k <- parms$k # m3/d
  ko <- parms$ko # cm/d
  #ka <- parms$ka ka = 0.089
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- c["Cw"]*(-k)
  # dmfdt:
  r[2] <- ko*Af*c["Cw"]/1000/L - ko*Af*c["mSPME"]/(Vf*Kf*L*1000) # Cw = [ng/L], mf = [ng/cm]
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{cinit <- c(Cw = 100000, mSPME = 0)
t <- c(1:40)
# Placeholder values of key parameters
parms <- list(k = 0.000, ko = 100)
out <- ode(y = cinit, times = t, func = rtm.SPME, parms = parms)
out <- data.frame(out)}

{new.out <- out %>%
  gather(variable, value, -time)
new.out <- within(new.out, variable <- factor(variable,
                                                levels = c('Cw', 'mSPME')))}

# Plot
ggplot(new.out, aes(x = time, y = value, color = variable)) +
  facet_wrap(vars(variable), scales = "free_y", nrow =  2) +
  geom_line(linewidth = 2) +
  theme_classic() +
  labs(x = 'time (day)', y = 'Concentration')
