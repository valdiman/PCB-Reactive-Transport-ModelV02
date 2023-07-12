# Code to model PCB4 in laboratory experiments
# using sediment from Altavista, VI. Passive measurements
# of PCB4 in the water and the air phases are predicted and
# linked to the water and air concentrations from the passive
# samplers.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("minpack.lm")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(minpack.lm) # least squares fit using levenberg-marquart algorithm
}

# Read data ---------------------------------------------------------------
exp.data.0 <- read.csv("Data/PCBDataV02.csv")

# Organize data -----------------------------------------------------------
# Remove lost sample(s), NA
exp.data <- exp.data.0[!is.na(exp.data.0$PCB4), ]

# spme = SPME fiber sampler
# puf = PUF sampler

# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
{exp.mspme.control <- exp.data %>%
  filter(treatment == "Ctrl" & sampler == "SPME") %>%
  group_by(time) %>%
  mutate(mean(PCB4/length)) %>%
  distinct(`mean(PCB4/length)`) %>%
  rename("PCB 4 mass in SPME (ng/cm); Control" = `mean(PCB4/length)`)

exp.mspme.LB400 <- exp.data %>%
  filter(treatment == "LB400" & sampler == "SPME") %>%
  group_by(time) %>%
  mutate(mean(PCB4/length)) %>%
  distinct(`mean(PCB4/length)`) %>%
  rename("PCB 4 mass in SPME (ng/cm); LB400" = `mean(PCB4/length)`) 

exp.mpuf.control <- exp.data %>%
  filter(treatment == "Ctrl" & sampler == "PUF") %>%
  group_by(time) %>%
  mutate(mean(PCB4)) %>%
  distinct(`mean(PCB4)`) %>%
  rename("PCB 4 mass in PUF (ng); Control" = `mean(PCB4)`)

exp.mpuf.LB400 <- exp.data %>%
  filter(treatment == "LB400" & sampler == "PUF") %>%
  group_by(time) %>%
  mutate(mean(PCB4)) %>%
  distinct(`mean(PCB4)`) %>%
  rename("PCB 4 mass in PUF (ng); LB400" = `mean(PCB4)`)

pcb.4 <- left_join(exp.mspme.control, exp.mpuf.control)  %>% 
  left_join(exp.mspme.LB400) %>% 
  left_join(exp.mpuf.LB400) %>% 
  mutate(time = recode(time, `1` = 16, `2` = 35, `3` = 75)) %>%
  ungroup() %>%
  add_row(time = 0, "PCB 4 mass in SPME (ng/cm); Control" = 0,
          "PCB 4 mass in PUF (ng); Control" = 0, 
          "PCB 4 mass in SPME (ng/cm); LB400" = 0,
          "PCB 4 mass in PUF (ng); LB400" = 0, .before = 1)}

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
  Va <- 125 # cm3 headspace volumne
  Vp <- 1 # cm3 (guest)
  Aaw <- 20 # cm2 
  Aws <- 30 # cm2
  
  # Congener-specific constants
  Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
  dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
  Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
  dUow <- -21338.96 # internal energy for the transfer of octanol-water for PCB 4 (J/mol)
  Koa <- 10^(6.521554861) # PCB 4 octanol-air equilibrium partition coefficient
  
  # PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366*log10(Koa) - 3.1774)# m3/g PCB 4-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 # g/m3 density of PUF
  # ro <- 0.0005 # m3/d sampling rate
  
  # SPME fiber constants
  Af <- 0.138 # cm2 SPME area
  Vf <- 0.000000069 # L/cm SPME volume/area
  L <- 30 # cm SPME length average
  Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
  # ko <- 70 # cm/d PCB 4 mass transfer coefficient to SPME
  
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
  
  # flux constant passed through a list called parms
  ro <- parms$ro # m3/d
  ko <- parms$ko # cm/d
  
  # derivatives dx/dt are computed below
  # C = ng/L
  r <- rep(0,length(c))
  # dCwdt:
  r[1] <- kaw.o*Aaw/Vw*(c["Ca"]/(Kaw.t) - c["Cw"]) + (D.pcb.water/bl)*Aws*60*60*24/Vw*(c["Cpw"] - c["Cw"]) + ko*Af/Vw*(c["mf"]/(Vf*Kf) - c["Cw"]) # 864 to change second to days and um to m
  # dCpwdt:
  r[2] <- (D.pcb.water/bl)*Aws*60*60*24/Vp*(c["Cw"]- c["Cpw"])
  # dCadt:
  r[3] <- kaw.o*Aaw/Va*(c["Cw"] - c["Ca"]/Kaw.t)*100^3 + ro/Va*(c["mpuf"]/(Vpuf*d*Kpuf*1000) - c["Ca"])
  # dmfdt:
  r[4] <- ko*Af*c["Cw"]/1000/L - ko*Af*c["mf"]/(Vf*L*Kf*1000) # Cw = [ng/L], mf = [ng/cm]
  # dmpufdt:
  r[5] <- ro*c["Ca"]*1000 - ro*(c["mpuf"]/(Vpuf*d))/(Kpuf) #  Ca = [ng/L], mpuf = [ng]
  
  # The computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

# Initial conditions and run function
{
  cinit <- c(Cw = 0, mf = 0, Ca = 0, mpuf = 0, Cpw = 242)
  t.1 <- pcb.4$time
  # Placeholder values of key parameters
  parms <- list(ro = 0.0056, ko = 5) # Input reasonable estimate of ko and ro (placeholder values)
  out.1 <- ode(y = cinit, times = t.1, func = rtm.PCB4, parms = parms)
  head(out.1)
}

# Fitting function --------------------------------------------------------
# Sums of squares function for parameter fitting 
ssq = function(parms){
  
  # initial concentration
  cinit <- c(Cw = 0, mf = 0, Ca = 0, mpuf = 0)
  
  # Time points for which values of mf and mpuf are reported
  # Include the points where data are available
  t <- c(seq(0, 75, 1), pcb.4$time)
  t <- sort(unique(t))
  
  # Parameters from the parameter estimation routine.
  # Values are forced to be positive to avoid errors in the solution
  ro <- sqrt((parms[1])^2) #
  ko <- sqrt((parms[2])^2) # 
  
  # Solve ODE for a given set of parameters
  out2 <- ode(y = cinit, times = t, func = rtm.PCB4,
              parms = list(ro = sqrt((ro)^2), ko = sqrt((ko)^2)))
  
  # Filter data that contains time points where data are available
  out2.pcb.4 <- data.frame(out2)
  out2.pcb.4 <- out2.pcb.4[out2.pcb.4$time %in% pcb.4$time, ]
  
  # Evaluate predicted vs experimental residual
  # spme and puf only 
  out2.pcb.4 <- subset(out2.pcb.4, select = -c(Cw, Ca))
  pred.pcb.4 <- melt(out2.pcb.4, id.var = "time", variable.name = "sampler",
                      value.name = "mass")
  # spme and puf only
  pcb.4.2 <- subset(pcb.4, select = c(time, `PCB 4 mass in SPME (ng/cm); Control`,
                                        `PCB 4 mass in PUF (ng); Control`))
  exp.pcb.4 <- melt(pcb.4.2,id.var = "time", variable.name = "sampler",
                     value.name = "mass")
  ssqres <- pred.pcb.4$mass - exp.pcb.4$mass
  
  # return predicted vs experimental residual
  return(ssqres)
}

# parameter fitting using levenberg marquart algorithm
parms <- c(ro = 0.000077, ko = 6.012) # Input second guess for parameters
# fitting
fitval <- nls.lm(par = parms, fn = ssq)
summary(fitval)
# Estimated parameter
parest <- as.list(coef(fitval))
parms <- list(ro = parest$ro, ko = parest$ko)

# Simulated predicted profile at estimated parameter values
{cinit <- c(Cw = 0, mf = 0, Ca = 0, mpuf = 0)
t.2 <- c(seq(0, 75, 1))
out3 <- ode(y = cinit, times = t.2, func = rtm.PCB4, parms = parest)
out.plot <- data.frame(out3)}

# Plotting ----------------------------------------------------------------
# (1) Plot of predicted vs experimental control data

names(out.plot) <- c("time", "PCB 4 predicted Cw", "PCB 4 predicted mass fiber",
                     "PCB 4 predicted Ca", "PCB 4 predicted mass PUF")
# Predicted
# spme
out.plot.mspme <- subset(out.plot,
                      select = c(time, `PCB 4 predicted mass fiber`))
# puf
out.plot.mpuf <- subset(out.plot,
                        select = c(time, `PCB 4 predicted mass PUF`))
pred.mspme <- melt(out.plot.mspme, id.var = c("time"),
                variable.name = "compartment", value.name = "mass")
pred.mpuf <- melt(out.plot.mpuf, id.var = c("time"),
                  variable.name = "compartment", value.name = "mass")
# Experimental
# spme
pcb.4.mspme <- subset(pcb.4, select = c(time, `PCB 4 mass in SPME (ng/cm); Control`))
# puf
pcb.4.mpuf <- subset(pcb.4, select = c(time, `PCB 4 mass in PUF (ng); Control`))
exp.mspme <- melt(pcb.4.mspme, id.var = c("time"), variable.name = "compartment",
               value.name = "mass")
exp.mpuf <- melt(pcb.4.mpuf, id.var = c("time"), variable.name = "compartment",
                 value.name = "mass")  
names(pcb.4.mspme) <- c("time", "PCB 4 measured mass fiber")
names(pcb.4.mpuf) <- c("time", "PCB 4 measured mass PUF")

mspme <- ggplot(data = pred.mspme, aes(x = time, y = mass)) +
  geom_line(colour = 2) +
  geom_point(data = exp.mspme, aes(x = time, y = mass), color = 2) +
  theme_bw() +
  theme(aspect.ratio = 6/10) +
  ggtitle("Control - SPME Fiber Results for PCB 4") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 4 fiber mass accumulated  (ng/cm)"))) +
  ylim(0, 0.25) +
  xlim(0, 80)

print(mspme)

mpuf <- ggplot(data = pred.mpuf, aes(x = time, y = mass)) +
  geom_line(colour = 4) +
  geom_point(data = exp.mpuf, aes(x = time, y = mass), color = 4) +
  theme_bw() +
  theme(aspect.ratio = 6/10) +
  ggtitle("Control - PUF Results for PCB 4") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 4 PUF mass accumulated (ng/PUF)"))) +
  ylim(0, 600) +
  xlim(0, 80)

print(mpuf)
