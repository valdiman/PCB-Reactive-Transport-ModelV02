# Load libraries
library(dplyr) # organize data
library(tidyverse)
library(ggplot2) # plotting
library(deSolve) # solving differential equations

# Reactive transport function ---------------------------------------------

rtm.PCB17 = function(t, c, parms){
  
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
  Kaw <- 0.015256157 # PCB 17 dimensionless Henry's law constant @ 25 C
  dUaw <- 52590.22 # internal energy for the transfer of air-water for PCB 4 (J/mol)
  Kow <- 10^(5.25) # PCB 17 octanol-water equilibrium partition coefficient
  Koa <- 10^(7.066554861) # PCB 17 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366*log10(Koa) - 3.1774)# m3/g PCB 4-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  
  # SPME fiber constants
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  
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
  Ct <- 307.3052312 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  logKoc <- 0.94*log10(Kow) + 0.42 # koc calculation
  Kd <- foc*10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/Kd*1000*(1) # [ng/L]
  
  # Biotransformation rate
  kb <- 0 # 0.2843 1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.163562
  
  # flux constant passed through a list called parms
  ro <- parms$ro # m3/d
  ko <- parms$ko # cm/d
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- kaw.o*Aaw/Vw*(c["Ca"]/(Kaw.t) - c["Cw"]) + D.pcb.water*Aws*60*60*24/bl/Vw*(Cpw - c["Cw"]) - kb*c["Cw"]
  # dmfdt:
  r[2] <- ko*Af*c["Cw"]/1000/L - ko*Af*c["mSPME"]/(Vf*Kf*L*1000) # Cw = [ng/L], mf = [ng/cm]
  # dCadt:
  r[3] <- kaw.o*Aaw/Va*(c["Cw"] - c["Ca"]/Kaw.t)
  # dmpufdt:
  r[4] <- ro*c["Ca"]*1000 - ro*(c["mPUF"]/(Vpuf*d))/(Kpuf) #  Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{cinit <- c(Cw = 0, mSPME = 0, Ca = 0, mPUF = 0)
t <- c(1:80)
# Placeholder values of key parameters
parms <- list(ro = 0.0045, ko = 100) # Input reasonable estimate of ko and ro (placeholder values)
out.PCB17 <- ode(y = cinit, times = t, func = rtm.PCB17, parms = parms) %>%
as.data.frame() -> out.PCB17
}

# Estimate % depletion from Cw
{L <- 30 # cm SPME length average
Vw <- 100 # cm3 water volume
out.PCB17$Depletion <- (out.PCB17$mSPME*L)/(out.PCB17$Cw*Vw/1000)*100}

# Estimate % depletion from total
{Ct <- 307.3052312 # ng/g PCB 17 sediment concentration
  ms <- 10 # g
  out.PCB17$DepletionT <- (out.PCB17$mSPME*L)/(Ct*ms)*100}

{new.out<- out.PCB17 %>%
  gather(variable, value, -time)
  new.out <- within(new.out,
                    variable <- factor(variable,
                                       levels = c('Cw', 'mSPME', 'Ca',
                                                  'mPUF', 'Depletion',
                                                  'DepletionT')))
}

# Plot
ggplot(new.out, aes(x = time, y = value, color = variable)) +
  facet_wrap(vars(variable), scales = "free_y", ncol = 2) +
  geom_line(size = 2) +
  theme_classic() +
  labs(x = 'time (day)', y = 'Concentration')

# Add experimental data
# Read data ---------------------------------------------------------------
exp.data <- read.csv("PCBDataV02.csv")

# Organize SPME data -----------------------------------------------------------
# spme = SPME fiber sampler [ng/cm]
# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
i <- "PCB17"
{exp.mspme <- exp.data %>%
    mutate(exp.data[all_of(i)]/length) %>%
    filter(sampler == "SPME") %>%
    group_by(time, treatment) %>%
    select(all_of(i)) %>%
    rename(!!paste(i, ".SPME", sep ="") := `i`) %>%
    mutate(time = recode(time, `1` = 16, `2` = 35, `3` = 75))
  
  exp.mspme.ctrl <- exp.mspme[1:9, 1:3] # experimental control control
  colnames(exp.mspme.ctrl) <- c('time', 'treatment', 'mSPME')
  exp.mspme.ctrl <- data.frame(exp.mspme.ctrl)
  exp.mspme.lb400 <- exp.mspme[10:18, 1:3] # experimental values LB400
  colnames(exp.mspme.lb400) <- c('time', 'treatment', 'mSPME')
  exp.mspme.lb400 <- data.frame(exp.mspme.lb400)
}

{out.mspme <- out.PCB17[, c(1,3)]
  out.mspme$treatment <- c('pred')
  out.mspme <- out.mspme %>% relocate(treatment, .before = mSPME) # predicted values
}

# mSPME plot
ggplot(NULL, aes(x = time, y = mSPME, color = treatment)) +
  geom_line(data = out.mspme) +
  geom_point(data = exp.mspme.ctrl, color = "blue") +
  geom_point(data = exp.mspme.lb400, color = "red") +
  theme_bw() +
  theme(aspect.ratio = 3/3) +
  labs(x = 'time (day)', y = paste(i, " SPME (ng/cm)", sep ="")) +
  scale_color_manual(values = c("ctrl"="blue", "lb400"="red",
                                "pred" = "black"))

# Organize PUF data -----------------------------------------------------------
# puf = PUF sampler [ng]
# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
{exp.mpuf <- exp.data %>%
  mutate(exp.data[all_of(i)]) %>%
  filter(sampler == "PUF") %>%
  group_by(time, treatment) %>%
  select(all_of(i)) %>%
  rename(!!paste(i, ".PUF", sep ="") := `i`) %>%
  mutate(time = recode(time, `1` = 16, `2` = 35, `3` = 75))

exp.mpuf.ctrl <- exp.mpuf[1:9, 1:3] # experimental control values
colnames(exp.mpuf.ctrl) <- c('time', 'treatment', 'mPUF')
exp.mpuf.ctrl <- data.frame(exp.mpuf.ctrl)
exp.mpuf.lb400 <- exp.mpuf[10:18, 1:3] # experimental values LB400
colnames(exp.mpuf.lb400) <- c('time', 'treatment', 'mPUF')
exp.mpuf.lb400 <- data.frame(exp.mpuf.lb400)

out.mpuf <- out.PCB17[, c(1,5)]
out.mpuf$treatment <- c('pred')
out.mpuf <- out.mpuf %>% relocate(treatment, .before = mPUF) # predicted values
}

# mPUF plot
ggplot(NULL, aes(x = time, y = mPUF, color = treatment)) +
  geom_line(data = out.mpuf) +
  geom_point(data = exp.mpuf.ctrl, color = "blue") +
  geom_point(data = exp.mpuf.lb400, color = "red") +
  theme_bw() +
  theme(aspect.ratio = 3/3) +
  labs(x = 'time (day)', y = paste(i, " PUF (ng)", sep ="")) +
  scale_color_manual(values = c("ctrl"="blue", "lb400"="red",
                                "pred" = "black"))

