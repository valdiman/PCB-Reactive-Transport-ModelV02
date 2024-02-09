
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("deSolve")

# Load libraries
{
  library(readxl) # read excel data
  library(dplyr) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
}

# Experimental data -------------------------------------------------------
# Read data from Excel sheets
data.1 <- read_excel("Data/SPMECalibration.xlsx", sheet = "day5_14")
data.2 <- read_excel("Data/SPMECalibration.xlsx", sheet = "day42")
data.3 <- read_excel("Data/SPMECalibration.xlsx", sheet = "day85")

# Calibration data --------------------------------------------------------
# Select individual congeners from datasets
pcbi <- "PCB4"

# Extract relevant columns from each dataset
{
  d.1.pcbi <- data.1[, c("sample", "treatment", "replicate",
                         "time", "length", pcbi)]
  d.2.pcbi <- data.2[, c("sample", "treatment", "replicate",
                         "time", "length", pcbi)]
  d.3.pcbi <- data.3[, c("sample", "treatment", "replicate",
                         "time", "length", pcbi)]
}

# Combine data frames into one
pcb.cali <- rbind(d.1.pcbi, d.2.pcbi, d.3.pcbi)

# Modify replicate variable to group specific levels together
pcb.cali$replicate_grouped <- ifelse(pcb.cali$replicate %in% c("r.3.1", "r.3.2", "r.3.3"), "r.3", pcb.cali$replicate)

# Recode the "r.3" levels in replicate_grouped to "r.3.n"
pcb.cali$replicate_grouped <- recode(pcb.cali$replicate_grouped,
                                     "r.3" = "r.3.n")

# Function of differential equations --------------------------------------
system_equations <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Unpack the state variables
    Cw <- state[1]  # Water concentration
    mf <- state[2]  # Additional variable
    
    # Constants
    Af <- 0.138 # cm2 SPME area
    Vf <- 0.000000069 # L/cm SPME volume/area
    L <- 30 # cm SPME length average
    Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
    Kf <- 10^(1.06*log10(Kow) - 1.16) # PCB 4-SPME equilibrium partition coefficient
     
    # Unpack the parameters
    kb <- rate_b  # Bio rate constant [1/d]
    kf <- mtcf.w  # MTC from fiber/water [cm/d]
    
    # Differential equations
    dCwdt <- - kb * Cw
    dmfdt <- (kf * Af / L / 1000) * (Cw - mf / (Kf * Vf))
    
    # Return the derivatives
    list(c(dCwdt, dmfdt))
  })
}

# Set initial conditions and parameters -----------------------------------
# Estimating Cpw (PCB 4 concentration in sediment porewater)
Ct <- 630.2023*(4.5) # ng/g PCB 4 sediment concentration
foc <- 0.03 # organic carbon % in sediment
Kow <- 10^(4.65) # PCB 4 octanol-water equilibrium partition coefficient
logKoc <- 0.94*log10(Kow) + 0.42 # koc calculation
Kd <- foc*10^(logKoc) # L/kg sediment-water equilibrium partition coefficient
Cpw <- Ct/Kd*1000 # [ng/L]
initial_conditions <- c(Cw = Cpw, mf = 0)
rate_b <- 0.03 # Bio rate constant
mtcf.w <- 80  # MTC from fiber/water

# Time vector
times <- seq(0, 100, by = 1)  # Time from 0 to 80 with intervals of 0.1

# Define parameters list
parameters <- list(rate_b = rate_b, mtcf.w = mtcf.w)

# Solve the system of differential equations ------------------------------
solution <- ode(y = initial_conditions, times = times,
                func = system_equations, parms = parameters)

# Plot --------------------------------------------------------------------
# Modeling Cw and mf
plot(solution, type = "l", xlab = "Time", ylab = "Concentration")

# Get max y-axis from experimental data
y_max <- max(pcb.cali$PCB4/pcb.cali$length)

# Plot the modeling results with specified y-axis limits
plot(solution[, "time"], solution[, "mf"], type = "l",
     xlab = "Time", ylab = "mf", ylim = c(0, y_max))

# Add experimental data points
points(pcb.cali$time, pcb.cali$PCB4/pcb.cali$length,
       col = "red", pch = 16)
