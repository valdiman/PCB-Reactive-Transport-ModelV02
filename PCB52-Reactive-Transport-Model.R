


# Install packages
install.packages("ggplot2")
install.packages("deSolve")
install.packages("minpack.lm")

# Load libraries
library(ggplot2) # plotting
library(deSolve) # solving differential equations
library(minpack.lm) # least squares fit using levenberg-marquart algorithm

# Read data.xlsx
# Mean values
pcb.52 <- read_excel("dataV03.xlsx", sheet = "pcb52",
                     col_names = TRUE, col_types = NULL)

# Reactive transport function ---------------------------------------------

rtm.PCB52 = function(t, c, parms){
  
  #Experimental conditions
  R <- 8.3144 #J/(mol K) molar gas constant 
  MH2O <- 18.0152 # g/mol water molecular weight
  Mco2 <- 44.0094 # g/mol CO2 molecular weight
  Tst <- 25 #C air temperature
  Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
  Tw <- 20 #C water temperature
  Tw.1 <- 273.15 + Tw
  
  ## bioreactor parameters
  Vpw <- 25/1000000 #m3 porewater volume
  Vw <- 100/1000000 #m3 water volume
  Va <- 125/1000000 #m3 headspace volumne
  Aaw <- 30 #cm2 for a ~3 cm radius air-water area
  
  #congener-specific constants
  Ka.w <-  0.0130452 # PCB52 dimensionless Henry's law constant @ 25 C
  dUow <- 55517.96 # internal energy for the transfer of octanol-water for PCB 52 (J/mol)
  logKoa <-  8.351339075 # PCB52 octanol-air equilibrium partition coefficient
  logKow <- 5.84 # PCB52 octanol-water equilibrium partition coefficient
  MW.pcb <- 291.976 # g/mol PCB 52 molecular weight
  
  #PUF constants 
  Vpuf <- 0.000029 # m3 volume of PUF
  Kpuf <- 10^(0.6366*logKoa-3.1774)# m3/g PCB 52-PUF equilibrium partition coefficient
  d <- 0.0213*100^3 #g/m3 density of PUF
  ro <- 0.0045 #m3/d sampling rate
  
  #SPME fiber constants
  Af <- 0.138 #cm2 SPME area
  Vf <- 0.000000069 #l/cm SPME volume/area
  L <- 20 #cm SPME length
  Kf <- 10^(1.06*logKow-1.16) # PCB 52-SPME equilibrium partition coefficient
  ko <- 70 #cm/d PCB 52 mass transfer coefficient to SPME
  
  #partitioning constants
  M <- 0.1 #kg/L solid-water ratio
  foc <- 0.03 # organic carbon % in particles
  K <- foc*(10^(0.94*logKow+0.42)) #L/kg sediment-water equilibrium partition coefficient
  D.water.air <- 0.2743615 # cm2/s water's diffusion coefficient in the gas phase @ Tair = 25 C, patm = 1013.25 mbars 
  D.co2.w <- 1.67606E-05 # cm2/s CO2's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars 
  D.pcb.air <- D.water.air*(MW.pcb/MH2O)^(-0.5) # cm2/s PCB 52's diffusion coefficient in the gas phase 
  D.pcb.water <- D.co2.w*(MW.pcb/Mco2)^(-0.5) # cm2/s PCB 52's diffusion coefficient in water @ Tair = 25 C, patm = 1013.25 mbars
  v.H2O <- 0.010072884	# cm2/s kinematic viscosity of water @ Tair = 25
  V.water.air <- 0.003 # m/s water's velocity of air-side mass transfer without ventilation
  V.co2.w <- 0.041 # m/s mass transfer coefficient of CO2 in water side without ventilation
  SC.pcb.w <- v.H2O/D.pcb.water # Schmidt number PCB 52
  
  # kaw calculations
  # i) Ka.w.t
  Ka.w.t <- Ka.w*exp(-dUow/R*(1/Tw.1-1/Tst.1))*Tst.1/Tw.1 # Ka.w corrected by water and air temps during experiment
  
  # ii) Kaw.a
  Kaw.a <- (D.pcb.air/D.water.air)^(0.67)*V.water.air # m/s air-side mass transfer coefficient
  
  # iii) Kaw.w
  Kaw.w <- V.co2.w*(SC.pcb.w/600)^(-0.5) # m/s water-side mass transfer coefficient for PCB 52. 600 is the Schmidt number of CO2 at 298 K
  
  # iv) kaw
  kaw <- (1/(Kaw.a*Ka.w.t) + (1/Kaw.w))^-1 # m/s overall air-water mass transfer coefficient for PCB 52
  kaw <- kaw*100*60*60*24 #cm/d overall air-water mass transfer coefficient for PCB 52
  
  # biotransformation rate
  kb <- 0 #1/d, value changes depending on experiment, i.e., control = 0, treatments LB400 = 0.130728499, LB400+Sap = 0.13325936
  
  # flux constant passed through a list called parms
  ka <- parms$ka #1/d
  kd <- parms$kd #1/d
  
  # derivatives dx/dt are computed below
  r <- rep(0,length(c))
  #r[1] <- c["Ca"]*Aaw*kaw/(Ka.w.t*Vw*10^6) + c["Cw"]*(kd*M*K*Vpw/Vw - ka - Aaw*kaw/(Vw*10^6)) -kb*c["Cw"] #dCwdt
  r[2] <- ko*Af*c["Cw"]/1000/L - ko*Af*c["mf"]/(Vf*L*Kf*1000) # dmfdt
  r[3] <- c["Cw"]*(kaw*Aaw/Va)/10^6 - c["Ca"]*(kaw*Aaw/Ka.w.t/Va)/10^6 # dCadt
  r[4] <- ro*c["Ca"]*1000 - ro*(c["mpuf"]/(Vpuf*d))/(Kpuf) # dmpufdt
  
  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
}

#Predicted initial PCB 52 concentration for a given parameter set
# Estimating Cwi (initial PCB 52 concentration in the water)
Ct <- 321.4900673 #ng/g PCB 52 sediment concentration
logKow <- 5.84 # PCB52 octanol-water equilibrium partition coefficient
foc <- 0.03 # organic carbon % in particles
K <- foc*(10^(0.94*logKow+0.42)) #L/kg PCB 52 sediment-water equilibrium partition coefficient
ds <- 900 #g/L density
M <- 0.1 #kg/L solid-water ratio
Cwi <- Ct*ds/(1+M*K)
cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
t <- pcb.52$time
# Placeholder values of key parameters
parms <- list(ka = 0.089, kd = 0.000036) #Input reasonable estimate of ka and kd (placeholder values)
out1 <- ode(y = cinit, times = t, func = rtm.PCB52, parms = parms)
head(out1)

#Sums of squares function for parameter fitting 
ssq = function(parms){
  
  # initial concentration
  cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
  
  # time points for which values of mf and mpuf are reported
  # include the points where data are available
  t <- c(seq(0, 35, 1), pcb.52$time)
  t <- sort(unique(t))
  
  # parameters from the parameter estimation routine. Values are forced to be positive to avoid errors in the solution
  ka <- sqrt((parms[1])^2) #1/d
  kd <- sqrt((parms[2])^2) #1/d
  
  # solve ODE for a given set of parameters
  out2 <- ode(y = cinit, times = t, func = rtm.PCB52,
              parms = list(ka = sqrt((ka)^2), kd = sqrt((kd)^2)))
  
  # Filter data that contains time points where data are available
  out2.pcb.52 <- data.frame(out2)
  out2.pcb.52 <- out2.pcb.52[out2.pcb.52$time %in% pcb.52$time,]
  
  # Evaluate predicted vs experimental residual
  out2.pcb.52 <- subset(out2.pcb.52, select = -c(Cw, Ca)) # only mass of fiber and puf
  pred.pcb.52 <- melt(out2.pcb.52, id.var = "time", variable.name = "sampler",
                      value.name = "mass")
  pcb.52.2 <- subset(pcb.52, select = c(time, `PCB 52 mass in SPME (ng/cm); Control`, `PCB 52 mass in PUF (ng); Control`)) # only mass of fiber and puf
  exp.pcb.52 <- melt(pcb.52.2,id.var = "time", variable.name = "sampler",
                     value.name = "mass")
  ssqres <- pred.pcb.52$mass-exp.pcb.52$mass
  
  # return predicted vs experimental residual
  return(ssqres)
}

# parameter fitting using levenberg marquart algorithm
parms <- c(ka = 14, kd = 0.0232) #Input second guess for parameters
# fitting
fitval <- nls.lm(par = parms, fn = ssq)
summary(fitval)
# Estimated parameter
parest <- as.list(coef(fitval))
parms <- list(ka = parest$ka, kd = parest$kd)

# simulated predicted profile at estimated parameter values
cinit <- c(Cw = Cwi, mf = 0, Ca = 0, mpuf = 0)
t.2 <- c(seq(0, 40, 0.5))
out3 <- ode(y = cinit, times = t.2, func = rtm.PCB52, parms = parest)
out.plot <- data.frame(out3)

# plot of predicted vs experimental data
# Overlay predicted profile with experimental data

names(out.plot) <- c("time", "PCB 52 predicted Cw", "PCB 52 predicted mass fiber",
                     "PCB 52 predicted Ca", "PCB 52 predicted mass PUF")
out.plot.mf <- subset(out.plot,
                      select = c(time, `PCB 52 predicted mass fiber`)) # only mass of fiber
out.plot.mpuf <- subset(out.plot,
                        select = c(time, `PCB 52 predicted mass PUF`)) # only mass of puf
pred.mf <- melt(out.plot.mf, id.var = c("time"),
                variable.name = "compartment", value.name = "mass")
pred.mpuf <- melt(out.plot.mpuf, id.var = c("time"),
                  variable.name = "compartment", value.name = "mass")
pcb.52.mf <- subset(pcb.52, select = c(time, `PCB 52 mass in SPME (ng/cm); Control`)) # only mass of fiber
pcb.52.mpuf <- subset(pcb.52, select = c(time, `PCB 52 mass in PUF (ng); Control`)) # only mass of puf
exp.mf <- melt(pcb.52.mf, id.var = c("time"), variable.name = "compartment",
               value.name = "mass")
exp.mpuf <- melt(pcb.52.mpuf, id.var = c("time"), variable.name = "compartment",
                 value.name = "mass")  
names(pcb.52.mf) <- c("time", "PCB 52 measured mass fiber")
names(pcb.52.mpuf) <- c("time", "PCB 52 measured mass PUF")

mf <- ggplot(data = pred.mf, aes(x = time, y = mass)) +
  geom_line(colour = 2) +
  geom_point(data = exp.mf, aes(x = time, y = mass), color = 2) +
  theme_bw() +
  theme(aspect.ratio = 6/10) +
  ggtitle("Control - SPME Fiber Results for PCB 52") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 52 fiber mass accumulated  (ng/cm)"))) +
  ylim(0, 0.4) +
  xlim(0, 40)

print(mf)

mpuf <- ggplot(data = pred.mpuf, aes(x = time, y = mass)) +
  geom_line(colour = 4) +
  geom_point(data = exp.mpuf, aes(x = time, y = mass), color = 4) +
  theme_bw() +
  theme(aspect.ratio = 6/10) +
  ggtitle("Control - PUF Results for PCB 52") +
  xlab(expression(bold("Time (day)"))) +
  ylab(expression(bold("PCB 52 PUF mass accumulated (ng/PUF)"))) +
  ylim(0, 80) +
  xlim(0, 40)

print(mpuf)
