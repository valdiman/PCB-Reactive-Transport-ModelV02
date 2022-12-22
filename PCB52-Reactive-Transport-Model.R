


# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("minpack.lm")

# Load libraries
library(dplyr) # organize data
library(reshape2) # organize data
library(ggplot2) # plotting
library(deSolve) # solving differential equations
library(minpack.lm) # least squares fit using levenberg-marquart algorithm

# Read data ---------------------------------------------------------------
exp.data.0 <- read.csv("PCBDataV02.csv")

# Organize data -----------------------------------------------------------
# Remove lost sample(s), NA
exp.data <- exp.data.0[!is.na(exp.data.0$PCB52), ]

# spme = SPME fiber sampler
# puf = PUF sampler

# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
exp.mspme.control <- exp.data %>%
  filter(treatment == "Ctrl" & sampler == "SPME") %>%
  group_by(time) %>%
  mutate(mean(PCB52/length)) %>%
  distinct(`mean(PCB52/length)`) %>%
  rename("PCB 52 mass in SPME (ng/cm); Control" = `mean(PCB52/length)`)

exp.mspme.LB400 <- exp.data %>%
  filter(treatment == "LB400" & sampler == "SPME") %>%
  group_by(time) %>%
  mutate(mean(PCB52/length)) %>%
  distinct(`mean(PCB52/length)`) %>%
  rename("PCB 52 mass in SPME (ng/cm); LB400" = `mean(PCB52/length)`) 

exp.mpuf.control <- exp.data %>%
  filter(treatment == "Ctrl" & sampler == "PUF") %>%
  group_by(time) %>%
  mutate(mean(PCB52)) %>%
  distinct(`mean(PCB52)`) %>%
  rename("PCB 52 mass in PUF (ng); Control" = `mean(PCB52)`)

exp.mpuf.LB400 <- exp.data %>%
  filter(treatment == "LB400" & sampler == "PUF") %>%
  group_by(time) %>%
  mutate(mean(PCB52)) %>%
  distinct(`mean(PCB52)`) %>%
  rename("PCB 52 mass in PUF (ng); LB400" = `mean(PCB52)`)

pcb.52 <- left_join(exp.mspme.control, exp.mpuf.control)  %>% 
  left_join(exp.mspme.LB400) %>% 
  left_join(exp.mpuf.LB400) %>% 
  mutate(time = recode(time, `1` = 16, `2` = 35, `3` = 75)) %>%
  ungroup() %>%
  add_row(time = 0, "PCB 52 mass in SPME (ng/cm); Control" = 0,
          "PCB 52 mass in PUF (ng); Control" = 0, 
          "PCB 52 mass in SPME (ng/cm); LB400" = 0,
          "PCB 52 mass in PUF (ng); LB400" = 0, .before = 1)

# Reactive transport function ---------------------------------------------

rtm.PCB52 = function(t, c, parms){
  
  # Experimental conditions
  R <- 8.3144 # J/(mol K) molar gas constant 
  MH2O <- 18.0152 # g/mol water molecular weight
  Mco2 <- 44.0094 # g/mol CO2 molecular weight
  Tst <- 25 # C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 13 # C water temperature
  Tw.1 <- 273.15 + Tw
  
  # Bioreactor parameters
  Vw <- 100/1000000 # m3 water volume
  Va <- 125/1000000 # m3 headspace volumne
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Ka.w <-  0.0130452 # PCB52 dimensionless Henry's law constant @ 25 C
  dUow <- 55517.96 # internal energy for the transfer of octanol-water for PCB 52 (J/mol)
  logKoa <-  8.351339075 # PCB52 octanol-air equilibrium partition coefficient
  logKow <- 5.84 # PCB52 octanol-water equilibrium partition coefficient
  MW.pcb <- 291.976 # g/mol PCB 52 molecular weight
  
  # PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366*logKoa - 3.1774)# m3/g PCB 52-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  # ro <- 0.0005 # m3/d sampling rate
  
  # SPME fiber constants
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  Kf <- 10^(1.06*logKow - 1.16) # PCB 52-SPME equilibrium partition coefficient
  # ko <- 70 # cm/d PCB 52 mass transfer coefficient to SPME
  
  # Air & water physical conditions
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 52's diffusion coefficient in the gas phase 
  D.pcb.water <- D.co2.w*(MW.pcb/Mco2)^(-0.5) # cm2/s PCB 52's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.001 # m/s water's velocity of air-side mass transfer without ventilation (eq. 20-15)
  V.co2.w <- 9*10^-6 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 52
  bl <- 0.2 # cm boundary layer thickness
  
  # kaw calculations (air-water mass transfer coefficient)
  # i) Ka.w.t, ka.w corrected by water and air temps during experiment
  Ka.w.t <- Ka.w*exp(-dUow/R*(1/Tw.1-1/Tst.1))*Tst.1/Tw.1
  # ii) Kaw.a, air-side mass transfer coefficient
  Kaw.a <- V.water.air*(D.pcb.air/D.water.air)^(0.67) # [m/s]
  # iii) Kaw.w, water-side mass transfer coefficient for PCB 52. 600 is the Schmidt number of CO2 at 298 K
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # [m/s] 
  # iv) kaw, overall air-water mass transfer coefficient for PCB 52
  kaw.o <- (1/(Kaw.a*Ka.w.t) + (1/Kaw.w))^-1 # [m/s]
  # v) kaw, overall air-water mass transfer coefficient for PCB 52, units change
  kaw.o <- kaw.o*100*60*60*24 # [cm/d]

  # Estimating Cpw (PCB 52 concentration in sediment porewater)
  Ct <- 321.4900673 # ng/g PCB 52 sediment concentration
  logKow <- 5.84 # PCB52 octanol-water equilibrium partition coefficient
  foc <- 0.03 # organic carbon % in sediment
  K <- foc*(10^(0.94*logKow + 0.42)) # L/kg sediment-water equilibrium partition coefficient
  Cpw <- Ct/K*1000 # [ng/L]
  
  # Biotransformation rate
  kb <- 0 # 1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.130728499

  # flux constant passed through a list called parms
  ro <- parms$ro # m3/d
  ko <- parms$ko # cm/d
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- kaw.o*Aaw/(Vw*10^6)*(c["Ca"]/(Ka.w.t) - c["Cw"]) - kb*c["Cw"] + D.pcb.water*Aws*0.0864/(bl*Vw)*(Cpw - c["Cw"]) # 864 to change second to days and um to m
  # dmfdt:
  r[2] <- ko*Af*c["Cw"]/1000/L - ko*Af*c["mf"]/(Vf*L*Kf*1000) # Cw = [ng/L], mf = [ng/cm]
  # dCadt:
  r[3] <- c["Cw"]*(kaw.o*Aaw/Va)/10^6 - c["Ca"]*(kaw.o*Aaw/Ka.w.t/Va)/10^6
  # dmpufdt:
  r[4] <- ro*c["Ca"]*1000 - ro*(c["mpuf"]/(Vpuf*d))/(Kpuf) #  Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
cinit <- c(Cw = 0, mf = 0, Ca = 0, mpuf = 0)
t.1 <- pcb.52$time
# Placeholder values of key parameters
parms <- list(ro = 0.0043, ko = 70) # Input reasonable estimate of ko and ro (placeholder values)
out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB52, parms = parms)
head(out.1)

# Fitting function --------------------------------------------------------
# Sums of squares function for parameter fitting 
ssq = function(parms){
  
  # initial concentration
  cinit <- c(Cw = 0, mf = 0, Ca = 0, mpuf = 0)
  
  # Time points for which values of mf and mpuf are reported
  # Include the points where data are available
  t <- c(seq(0, 75, 1), pcb.52$time)
  t <- sort(unique(t))
  
  # Parameters from the parameter estimation routine.
  # Values are forced to be positive to avoid errors in the solution
  ro <- sqrt((parms[1])^2) #
  ko <- sqrt((parms[2])^2) # 
  
  # Solve ODE for a given set of parameters
  out2 <- ode(y = cinit, times = t, func = rtm.PCB52,
              parms = list(ro = sqrt((ro)^2), ko = sqrt((ko)^2)))
  
  # Filter data that contains time points where data are available
  out2.pcb.52 <- data.frame(out2)
  out2.pcb.52 <- out2.pcb.52[out2.pcb.52$time %in% pcb.52$time, ]
  
  # Evaluate predicted vs experimental residual
  # spme and puf only 
  out2.pcb.52 <- subset(out2.pcb.52, select = -c(Cw, Ca))
  pred.pcb.52 <- melt(out2.pcb.52, id.var = "time", variable.name = "sampler",
                      value.name = "mass")
  # spme and puf only
  pcb.52.2 <- subset(pcb.52, select = c(time, `PCB 52 mass in SPME (ng/cm); Control`,
                                        `PCB 52 mass in PUF (ng); Control`))
  exp.pcb.52 <- melt(pcb.52.2,id.var = "time", variable.name = "sampler",
                     value.name = "mass")
  ssqres <- pred.pcb.52$mass - exp.pcb.52$mass
  
  # return predicted vs experimental residual
  return(ssqres)
}

# parameter fitting using levenberg marquart algorithm
parms <- c(ro = 0.0075, ko = 20) # Input second guess for parameters
# fitting
fitval <- nls.lm(par = parms, fn = ssq)
summary(fitval)
# Estimated parameter
parest <- as.list(coef(fitval))
parms <- list(ro = parest$ro, ko = parest$ko)

# Simulated predicted profile at estimated parameter values
cinit <- c(Cw = 0, mf = 0, Ca = 0, mpuf = 0)
t.2 <- c(seq(0, 75, 1))
out3 <- ode(y = cinit, times = t.2, func = rtm.PCB52, parms = parest)
out.plot <- data.frame(out3)

# Plotting ----------------------------------------------------------------
# (1) Plot of predicted vs experimental control data

names(out.plot) <- c("time", "PCB 52 predicted Cw", "PCB 52 predicted mass fiber",
                     "PCB 52 predicted Ca", "PCB 52 predicted mass PUF")
# Predicted
# spme
out.plot.mspme <- subset(out.plot,
                      select = c(time, `PCB 52 predicted mass fiber`))
# puf
out.plot.mpuf <- subset(out.plot,
                        select = c(time, `PCB 52 predicted mass PUF`))
pred.mspme <- melt(out.plot.mspme, id.var = c("time"),
                variable.name = "compartment", value.name = "mass")
pred.mpuf <- melt(out.plot.mpuf, id.var = c("time"),
                  variable.name = "compartment", value.name = "mass")
# Experimental
# spme
pcb.52.mspme <- subset(pcb.52, select = c(time, `PCB 52 mass in SPME (ng/cm); Control`))
# puf
pcb.52.mpuf <- subset(pcb.52, select = c(time, `PCB 52 mass in PUF (ng); Control`))
exp.mspme <- melt(pcb.52.mspme, id.var = c("time"), variable.name = "compartment",
               value.name = "mass")
exp.mpuf <- melt(pcb.52.mpuf, id.var = c("time"), variable.name = "compartment",
                 value.name = "mass")  
names(pcb.52.mspme) <- c("time", "PCB 52 measured mass fiber")
names(pcb.52.mpuf) <- c("time", "PCB 52 measured mass PUF")

mspme <- ggplot(data = pred.mspme, aes(x = time, y = mass)) +
  geom_line(colour = 2) +
  geom_point(data = exp.mspme, aes(x = time, y = mass), color = 2) +
  theme_bw() +
  theme(aspect.ratio = 6/10) +
  ggtitle("Control - SPME Fiber Results for PCB 52") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 52 fiber mass accumulated  (ng/cm)"))) +
  ylim(0, 0.25) +
  xlim(0, 80)

print(mspme)

mpuf <- ggplot(data = pred.mpuf, aes(x = time, y = mass)) +
  geom_line(colour = 4) +
  geom_point(data = exp.mpuf, aes(x = time, y = mass), color = 4) +
  theme_bw() +
  theme(aspect.ratio = 6/10) +
  ggtitle("Control - PUF Results for PCB 52") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 52 PUF mass accumulated (ng/PUF)"))) +
  ylim(0, 40) +
  xlim(0, 80)

print(mpuf)
