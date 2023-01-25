# Calculation of the fraction of PCBs in the
# various compartments, sediment, water, sampler, air

# Constant ----------------------------------------------------------------
ms <- 10 # g
foc <- 0.03 # organic carbon % in sediment kgoc/kg
Vw <- 100/1000 # L water volume
Va <- 125/1000 # L head-space volume
Vf <- 0.000000069 # L/cm SPME volume/length
L <- 30 # cm SPME length average
Tst <- 25 # C air temperature
Tst.1 <- 273.15 + Tst # air and standard temperature in K, 25 C
Tw <- 13 # C water temperature
Tw.1 <- 273.15 + Tw # water temperature in K
R <- 8.3144 # J/(mol K) molar gas constant
# PCB 4
Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient Lw/kgoctanol
Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient Lw/LSPME
logKoc <- 0.94*log10(Kow) + 0.42 # koc calculation Lw/kgoc
Koc <- 10^(logKoc) # Lw/kgoc

# Calculations ------------------------------------------------------------
# Kaw.t, ka.w corrected by water and air temps during experiment
Kaw.t <- Kaw*exp(-dUaw/R*(1/Tw.1-1/Tst.1))*Tw.1/Tst.1
# Non-depletion calculation
ND <- Kf*Vf*L/Vw
ND < 1
Vf*L*Kf < Koc*ms*foc/1000
Vf*L*Kf/(Koc*ms*foc/1000) # much less than 1

# Fractions @ equilibrium
# (i) Sediment, water
den1 <- 1 + ((Koc*ms*foc/1000)/Vw)
fw1 <- 1/den1
fsed1 <- ((Koc*ms*foc/1000)/Vw)/den1

# (ii) Sediment, water, SPME
den2 <- 1 + ((Kf*Vf*L)/Vw) + ((Koc*ms*foc/1000)/Vw)
fw2 <- 1/den2
fsed2 <- ((Koc*ms*foc/1000)/Vw)/den2
fSPME <- ((Kf*Vf*L)/Vw)/den2
# Mass balance
fSPME/fw2*100
fSPME/fsed2*100
fw2/fsed2*100

# (iii) Sediment, water, air
den3 <- 1 + ((Kf*Vf*L)/Vw) + ((Koc*ms*foc/1000)/Vw) + (Kaw.t*Va/Vw)
fw3 <- 1/den3
fsed3 <- ((Koc*ms*foc/1000)/Vw)/den3
fSPME2 <- ((Kf*Vf*L)/Vw)/den3
fa1 <- (Kaw.t*Va/Vw)/den3
fSPME2/fw3*100
