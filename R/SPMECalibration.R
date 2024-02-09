
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")

# Load libraries
{
  library(readxl) # read excel data
  library(dplyr) # organize data
  library(ggplot2) # plotting
}

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
    K <- 5000 # Partition coefficient
     
    # Unpack the parameters
    kb <- rate_b  # Bio rate constant [1/d]
    kf <- mtcf.w  # MTC from fiber/water [cm/d]
    
    # Differential equations
    dCwdt <- - kb * Cw
    dmfdt <- (kf * Af / L / 1000) * (Cw - mf / (K * Vf))
    
    # Return the derivatives
    list(c(dCwdt, dmfdt))
  })
}

# Set initial conditions and parameters -----------------------------------
initial_conditions <- c(Cw = 300, mf = 0)
rate_b <- 0.3 # Bio rate constant
mtcf.w <- 5  # MTC from fiber/water

# Time vector
times <- seq(0, 80, by = 1)  # Time from 0 to 10 with intervals of 0.1

# Define parameters list
parameters <- list(rate_b = rate_b, mtcf.w = mtcf.w)

# Solve the system of differential equations ------------------------------
solution <- ode(y = initial_conditions, times = times,
                func = system_equations, parms = parameters)

# Plot --------------------------------------------------------------------
plot(solution, type = "l", xlab = "Time", ylab = "Concentration")


