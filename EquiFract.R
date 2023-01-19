# Calculation of the fraction of PCBs in the
# various compartments, sediment, water, sampler, air

# Not depletion
Vw <- 100/1000 # L water volume
Va <- 125/1000 # L head-space volume
Vf <- 0.000000069 # L/cm SPME volume/area
L <- 30 # cm SPME length average
Kaw <- 0.01344142 # PCB 4 dimensionless Henry's law constant @ 25 C
dUaw <- 49662.48 # internal energy for the transfer of air-water for PCB 4 (J/mol)
Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient Lw/kgoctanol
Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient Lw/LSPME

ND <- Kf*Vf*L/Vw
ND < 1

# Fractions
Ct <- 630.2023 # ng/g PCB 4 sediment concentration
ms <- 10 # g
foc <- 0.03 # organic carbon % in sediment kgoc/kg
logKoc <- 0.94*log10(Kow) + 0.42 # koc calculation L/kgoc
Koc <- 10^(logKoc) # Lw/kgoc

# (i) Sediment, water
den1 <- (1 + 1*(Koc*ms*foc/1000)/Vw)
fw1 <- 1/den1
fsed1 <- 1*(Koc*ms*foc/1000)/Vw/den1

# (ii) Sediment, water, SPME
den2 <- (1 + 1*(Kf*Vf*L)/Vw + 1*(Koc*ms*foc/1000)/Vw)
fw2 <- 1/den2
fsed2 <- 1*(Koc*ms*foc/1000)/Vw/den2
fSPME <- 1*(Kf*Vf*L)/Vw/den2
# Mass balance
mt <- 100 # it doesn't matter the units
mw1 <- mt*fw1
mw2 <- mt*fw2
mw1-mw2
mSPME <- mt*fSPME


