# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")

# Load libraries
library(dplyr) # organize data
library(reshape2) # organize data
library(ggplot2) # plotting


# Read data ---------------------------------------------------------------
exp.data.0 <- read.csv("PCBDataV02.csv")

# Organize data -----------------------------------------------------------
# Remove lost sample(s), NA
exp.data <- exp.data.0
exp.data <- exp.data.0[!is.na(exp.data.0$PCB4), ]

# spme = SPME fiber sampler
# puf = PUF sampler

# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
exp.mspme.control <- exp.data %>%
  filter(treatment == "Ctrl" & sampler == "SPME") %>%
  group_by(time) %>%
  mutate(PCB4/length) %>%
  distinct(`PCB4/length`) %>%
  rename("PCB 4 mass in SPME (ng/cm); Control" = `PCB4/length`)

exp.mspme.LB400 <- exp.data %>%
  filter(treatment == "LB400" & sampler == "SPME") %>%
  group_by(time) %>%
  mutate(PCB4/length) %>%
  distinct(`PCB4/length`) %>%
  rename("PCB 4 mass in SPME (ng/cm); LB400" = `PCB4/length`) 

exp.mpuf.control <- exp.data %>%
  filter(treatment == "Ctrl" & sampler == "PUF") %>%
  group_by(time) %>%
  mutate(PCB52) %>%
  distinct(`PCB52`) %>%
  rename("PCB 4 mass in PUF (ng); Control" = `PCB52`)

exp.mpuf.LB400 <- exp.data %>%
  filter(treatment == "LB400" & sampler == "PUF") %>%
  group_by(time) %>%
  mutate(PCB52) %>%
  distinct(`PCB52`) %>%
  rename("PCB 4 mass in PUF (ng); LB400" = `PCB52`)

PCB.4 <- cbind(exp.mspme.control, exp.mspme.LB400$`PCB 4 mass in SPME (ng/cm); LB400`,
               exp.mpuf.control$`PCB 4 mass in PUF (ng); Control`, exp.mpuf.LB400$`PCB 4 mass in PUF (ng); LB400`)
pCB.4 <-   add_row(time = 0, "PCB 4 mass in SPME (ng/cm); Control" = 0,
                   "PCB 4 mass in PUF (ng); Control" = 0, 
                   "PCB 4 mass in SPME (ng/cm); LB400" = 0,
                   "PCB 4 mass in PUF (ng); LB400" = 0, .before = 1)

a <- inner_join(exp.mspme.control, exp.mspme.LB400, by = "time")


a <- left_join(exp.mspme.control %>% group_by('time') %>% mutate(id = row_number()),
          exp.mspme.LB400 %>% group_by('time') %>% mutate(id = row_number()))


a <- mutate(exp.mspme.control, id = row_number())

a <- left_join(exp.mspme.control, exp.mspme.LB400) %>%
mutate('time')

a <- merge(exp.mspme.control, exp.mspme.LB400, by = 'time')

  exp.mspme.control %>%
  full_join(exp.mspme.LB400)


pcb.4 <- left_join(exp.mspme.control, exp.mpuf.control)  %>% 
  left_join(exp.mspme.LB400) %>% 
  left_join(exp.mpuf.LB400) %>% 
  mutate(time = recode(time, `1` = 16, `2` = 35, `3` = 75)) %>%
  ungroup() %>%
  add_row(time = 0, "PCB 4 mass in SPME (ng/cm); Control" = 0,
          "PCB 4 mass in PUF (ng); Control" = 0, 
          "PCB 4 mass in SPME (ng/cm); LB400" = 0,
          "PCB 4 mass in PUF (ng); LB400" = 0, .before = 1)

a <- left_join(exp.mspme.control, exp.mspme.LB400)

