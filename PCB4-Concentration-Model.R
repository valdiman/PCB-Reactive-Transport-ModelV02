# Load libraries
library(dplyr) # organize data
library(tidyverse)
library(ggplot2) # plotting
library(deSolve) # solving differential equations

# Reactive transport function ---------------------------------------------

rtm.PCB4 = function(t, c, parms){
  
  # Experimental conditions
  MH2O <- 18.0152 # g/mol water molecular weight
  MCO2 <- 44.0094 # g/mol CO2 molecular weight
  MW.pcb <- 223.088 # g/mol PCB 4 molecular weight
  Tst <- 25 # C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 13 # C water temperature
  Tw.1 <- 273.15 + Tw # water temperature in K
  R <- 8.3144 # J/(mol K) molar gas constant 
  
  # Bioreactor parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 head-space volume
  Aaw <- 20 # cm2 air-water area
  Aws <- 30 # cm2 sediment-water area
  
  # Congener-specific constants
  Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
  dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <- -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 4's diffusion coefficient in the gas phase (eq. 18-45)
  D.pcb.water <- D.co2.w*(MW.pcb/MCO2)^(-0.5) # cm2/s PCB 4's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 4.1*10^-2 # m/s mass transfer coefficient of CO2 in water side without ventilation
  # V.co2.w.2 <- 9*10^(-4)/100 # m/s new book 19-20
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 4
  bl <- 0.21 # cm boundary layer thickness
  
  # kaw calculations (air-water mass transfer coefficient)
  # i) Ka.w.t, ka.w corrected by water and air temps during experiment
  Kaw.t <- Kaw*exp(-dUaw/R*(1/Tw.1-1/Tst.1))*Tw.1/Tst.1
  # ii) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # iii) Kaw.w, water-side mass transfer coefficient for PCB 4. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s]
  # Kaw.w.2 <- V.co2.w.2*(SC.pcb.w/600)^(-0.5) # [m/s]
  # iv) kaw, overall air-water mass transfer coefficient for PCB 4
  kaw.o <- (1/(Kaw.a*Kaw.t) + (1/Kaw.w))^-1 # [m/s]
  # v) kaw, overall air-water mass transfer coefficient for PCB 4, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]
  
  # Estimating Cpw (PCB 4 concentration in sediment porewater)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kow.t <- Kow*exp(-dUow/R*(1/Tw.1-1/Tst.1)) # temperature correction
  logKoc <- 0.94*log10(Kow.t) + 0.42 # koc calculation
  Kd <- foc*10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/Kd*1000 # [ng/L]
  
  # Biotransformation rate
  kb <- 0 #0.130728499 # 1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.130728499
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- kaw.o*Aaw/Vw*(c["Ca"]/(Kaw.t) - c["Cw"]) + D.pcb.water*Aws*60*60*24/bl/Vw*(Cpw - c["Cw"]) - kb*c["Cw"]
   # dCadt:
  r[2] <- kaw.o*Aaw/Va*(c["Cw"] - c["Ca"]/Kaw.t)
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{cinit <- c(Cw = 0, Ca = 0)
t.1 <- c(1:10)
# Placeholder values of key parameters
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4) %>%
as.data.frame() -> out.1
}

{new.out<- out.1 %>%
  gather(variable, value, -time)
new.out <- within(new.out, variable <- factor(variable,
                                              levels = c('Cw', 'Ca')))
}

# Plot
ggplot(new.out, aes(x = time, y = value, color = variable)) +
  facet_wrap(vars(variable), scales = "free_y", ncol = 2) +
  geom_line(size = 2) +
  theme_classic() +
  labs(x = 'time (day)', y = 'Concentration')
