# Load libraries
library(dplyr) # organize data
library(tidyverse)
library(ggplot2) # plotting
library(deSolve) # solving differential equations

# Reactive transport function ---------------------------------------------

rtm.SPME.1 = function(t, c, parms){
  
  Vw <- 100/1000 # L water volume
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  
  # Estimating Cpw (PCB 4 concentration in sediment porewater)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  logKoc <- 0.94*log10(Kow) + 0.42 # koc calculation
  Kd <- foc*10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/Kd*1000*(1) # [ng/L]
  
  # Constant passed through a list called parms
  k <- parms$k # m3/d
  ko <- parms$ko # cm/d
  #ka <- parms$ka ka = 0.089
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- c["Cw"]*k + 10
  # dmfdt:
  r[2] <- ko*Af*c["Cw"]/1000/L - ko*Af*c["mSPME"]/(Vf*Kf*L*1000) # Cw = [ng/L], mf = [ng/cm]
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{cinit <- c(Cw = 0, mSPME = 0)
t <- c(1:25)
# Placeholder values of key parameters
parms <- list(k = -0.2, ko = 10)
out.1 <- ode(y = cinit, times = t, func = rtm.SPME.1, parms = parms)
out.1 <- data.frame(out.1)}

{new.out.1 <- out.1 %>%
  gather(variable, value, -time)
new.out.1 <- within(new.out.1, variable <- factor(variable,
                                                levels = c('Cw', 'mSPME')))}


rtm.SPME.2 = function(t, c, parms){
  
  Vw <- 100/1000 # L water volume
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  
  # Estimating Cpw (PCB 4 concentration in sediment porewater)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  logKoc <- 0.94*log10(Kow) + 0.42 # koc calculation
  Kd <- foc*10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/Kd*1000*(1) # [ng/L]
  
  # Constant passed through a list called parms
  k <- parms$k # m3/d
  ko <- parms$ko # cm/d
  #ka <- parms$ka ka = 0.089
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- c["Cw"]*k + 10
  # dmfdt:
  r[2] <- - ko*Af*c["mSPME"]/(Vf*Kf*L*1000)
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{cinit <- c(Cw = out.1$Cw[25], mSPME = out.1$mSPME[25])
  t <- c(26:80)
  # Placeholder values of key parameters
  parms <- list(k = -0.2, ko = 10)
  out.2 <- ode(y = cinit, times = t, func = rtm.SPME.2, parms = parms)
  out.2 <- data.frame(out.2)}

{new.out.2 <- out.2 %>%
    gather(variable, value, -time)
  new.out.2 <- within(new.out.2, variable <- factor(variable,
                                                    levels = c('Cw', 'mSPME')))}
# Merge both results (new.out.1 and mew.out.2)
new.out.3 <- merge(x = new.out.1, y = new.out.2, all = TRUE)

# Plot
ggplot(new.out.3, aes(x = time, y = value, color = variable)) +
  facet_wrap(vars(variable), scales = "free_y", nrow =  2) +
  geom_line(linewidth = 2) +
  theme_classic() +
  labs(x = 'time (day)', y = 'Concentration')
