# Load libraries
library(dplyr) # organize data
library(tidyverse)
library(ggplot2) # plotting
library(deSolve) # solving differential equations

# Reactive transport function ---------------------------------------------

rtm.SPME.1 = function(t, c, parms){
  
  # Parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 head-space volume
  Aaw <- 20 # cm2 air-water area
  Aws <- 30 # cm2 sediment-water area
  Vpuf <- 0.000029 # m3 volume of PUF
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  d <- 0.0213*100^3 # g/m3 density of PUF
  
  # Estimating Cpw (PCB 4 concentration in sediment porewater)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kd <- foc*10^(4.791) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/Kd*1000*(1.2) # [ng/L]
  
  # Partitioning and mass transfer
  kaw.o <- 62.1094 # [cm/d]
  Kaw.t <- 0.0055684
  D.pcb.water <- 7.444305e-06 # [cm2/s]
  Kpuf <- 9.423708
  Kf <- 5874.894
  bl <- 0.21 # cm boundary layer thickness
  
  # Biotransformation rate
  kb <- 0 #0.6161 # 1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.178374771
  
  # flux constant passed through a list called parms
  ko <- parms$ko # cm/d
  ro <- parms$ro # m3/d
  
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
t <- c(1:23)
# Placeholder values of key parameters
parms <- list(ko = 10, ro = 0.005)
out.1 <- ode(y = cinit, times = t, func = rtm.SPME.1, parms = parms)
out.1 <- data.frame(out.1)}

{new.out.1 <- out.1 %>%
  gather(variable, value, -time)
new.out.1 <- within(new.out.1, variable <- factor(variable,
                                                levels = c('Cw', 'mSPME', 'Ca', 'mPUF')))}

rtm.SPME.2 = function(t, c, parms){
  
  # Parameters
  Vw <- 100 # cm3 water volume
  Va <- 125 # cm3 head-space volume
  Aaw <- 20 # cm2 air-water area
  Aws <- 30 # cm2 sediment-water area
  Vpuf <- 0.000029 # m3 volume of PUF
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  d <- 0.0213*100^3 # g/m3 density of PUF
  
  # Estimating Cpw (PCB 4 concentration in sediment porewater)
  Ct <- 630.2023 # ng/g PCB 4 sediment concentration
  foc <- 0.03 # organic carbon % in sediment
  Kd <- foc*10^(4.791) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/Kd*1000*(1.2) # [ng/L]
  
  # Partitioning and mass transfer
  kaw.o <- 62.1094 # [cm/d]
  Kaw.t <- 0.0055684
  D.pcb.water <- 7.444305e-06 # [cm2/s]
  Kpuf <- 9.423708
  Kf <- 5874.894
  bl <- 0.21 # cm boundary layer thickness
  
  # Biotransformation rate
  kb <- 0 # 0.6161 # 1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.178374771
  
  # flux constant passed through a list called parms
  ko <- parms$ko # cm/d
  ro <- parms$ro # m3/d
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- kaw.o*Aaw/Vw*(c["Ca"]/(Kaw.t) - c["Cw"]) + D.pcb.water*Aws*60*60*24/bl/Vw*(Cpw - c["Cw"]) - kb*c["Cw"]
  # dmfdt:
  r[2] <- - ko*Af*c["mSPME"]/(Vf*Kf*L*1000)
  # dCadt:
  r[3] <- kaw.o*Aaw/Va*(c["Cw"] - c["Ca"]/Kaw.t)
  # dmpufdt:
  r[4] <- ro*c["Ca"]*1000 - ro*(c["mPUF"]/(Vpuf*d))/(Kpuf) #  Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{cinit <- c(Cw = out.1$Cw[23], mSPME = out.1$mSPME[23], Ca = out.1$Ca[23], mPUF = out.1$mPUF[23]) # 90% of mSPME
  t <- c(24:80)
  # Placeholder values of key parameters
  parms <- list(ko = 10, ro = 0.005)
  out.2 <- ode(y = cinit, times = t, func = rtm.SPME.2, parms = parms)
  out.2 <- data.frame(out.2)}

{new.out.2 <- out.2 %>%
    gather(variable, value, -time)
  new.out.2 <- within(new.out.2, variable <- factor(variable,
                                                    levels = c('Cw', 'mSPME', 'Ca', 'mPUF')))}
# Merge both results (new.out.1 and mew.out.2)
new.out.3 <- merge(x = new.out.1, y = new.out.2, all = TRUE)

# Plot
ggplot(new.out.3, aes(x = time, y = value, color = variable)) +
  facet_wrap(vars(variable), scales = "free_y", nrow =  2) +
  geom_line(linewidth = 2) +
  theme_classic() +
  labs(x = 'time (day)', y = 'Concentration')

# Add experimental data
# Read data ---------------------------------------------------------------
exp.data <- read.csv("PCBDataV02.csv")

# Organize SPME data -----------------------------------------------------------
# spme = SPME fiber sampler [ng/cm]
# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
i <- "PCB4"
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

out.3 <- merge(out.1, out.2, all = TRUE)

{out.mspme <- out.3[, c(1,3)]
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

out.mpuf <- out.3[, c(1,5)]
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



