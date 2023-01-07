# Statistical analysis to review data, mostly between control
# and experiments, same time points.

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
exp.data <- read.csv("PCBDataV02.csv")

# Organize data -----------------------------------------------------------
# Remove lost sample(s), NA
# Not sure, it deoends on the t test
#exp.data <- exp.data.0[!(rowSums(exp.data.0[, c(8:180)],
#                                   na.rm = TRUE)==0),]

# spme = SPME fiber sampler [ng/cm]
# puf = PUF sampler [ng]

# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
# (i) Select congener, i
i <- "PCB4"
exp.mspme <- exp.data %>%
  mutate(exp.data[i]/length) %>%
  filter(sampler == "SPME") %>%
  group_by(time, treatment) %>%
  select(all_of(i)) %>%
  rename("PCB4.SPME" = `i`)
# (ii) Organize data for t test
# Select time, t
t <- "1"
exp.mspme.t.1 <- exp.mspme %>%
  filter(time == t)
exp.mspme.t.1 <- data.frame(Column1 = exp.mspme.t.1$`PCB4.SPME`[c(1:3)], 
           Column2 = exp.mspme.t.1$`PCB4.SPME`[c(4:6)])
exp.mspme.t.1 <- cbind(c(t, t, t), exp.mspme.t.1$Column1,
                       exp.mspme.t.1$Column2)
exp.mspme.t.1 <- as.numeric(exp.mspme.t.1)
colnames(exp.mspme.t.1) <- c("time", "PCB4.SPME[Control]",
                             "PCB4.SPME[LB400]")

# t-test
t.test(exp.mspme.t.1$`PCB4.SPME[Control]`,
       exp.mspme.t.1$`PCB4.SPME[LB400]`, var.equal = TRUE)


# LB400
exp.mspme.t.1.LB400 <- exp.mspme %>%
  filter(time == t & treatment == "LB400")
# Combine
pcb4.t.1 <- data.frame(cbind(exp.mspme.t.1.ctrl[1],
                   exp.mspme.t.1.ctrl[3],
                   exp.mspme.t.1.LB400[3]))
# Add column names
colnames(pcb4.t.1) <- c("time", "PCB4 SPME (ng/cm); Control",
                         "PCB4 SPME (ng/cm); LB400")
# t-test
t.test(pcb4.t.1$`PCB4.SPME; Control`,
       pcb4.t.1$`PCB4.SPME; LB400`, var.equal = TRUE)

# Time 2 and ctrl
exp.mspme.t.2.ctrl <- exp.mspme %>%
  filter(time == 2 & treatment == "Ctrl") %>%
  rename("Ctrl" = 'PCB4.SPME.(ng/cm)')
# Time 2 and LB400
exp.mspme.t.2.LB400 <- exp.mspme %>%
  filter(time == 2 & treatment == "LB400") %>%
  rename("LB400" = 'PCB4.SPME.(ng/cm)')
# Combine
pcb4.t.2 <- data.frame(cbind(exp.mspme.t.2.ctrl$time,
                             exp.mspme.t.2.ctrl$Ctrl,
                             exp.mspme.t.2.LB400$LB400))
# Add column names
colnames(pcb4.t.2) <- c("time", "PCB4 SPME (ng/cm); Control",
                        "PCB4 SPME (ng/cm); LB400")
# t-test
t.test(pcb4.t.2$`PCB4 SPME (ng/cm); Control`,
       pcb4.t.2$`PCB4 SPME (ng/cm); LB400`)

# Time 3 and ctrl
exp.mspme.t.3.ctrl <- exp.mspme %>%
  filter(time == 3 & treatment == "Ctrl") %>%
  rename("Ctrl" = 'PCB4.SPME.(ng/cm)')
# Time 3 and LB400
exp.mspme.t.3.LB400 <- exp.mspme %>%
  filter(time == 3 & treatment == "LB400") %>%
  rename("LB400" = 'PCB4.SPME.(ng/cm)')
# Combine
pcb4.t.3 <- data.frame(cbind(exp.mspme.t.3.ctrl$time,
                             exp.mspme.t.3.ctrl$Ctrl,
                             exp.mspme.t.3.LB400$LB400))
# Add column names
colnames(pcb4.t.3) <- c("time", "PCB4 SPME (ng/cm); Control",
                        "PCB4 SPME (ng/cm); LB400")
# t-test
t.test(pcb4.t.3$`PCB4 SPME (ng/cm); Control`,
       pcb4.t.3$`PCB4 SPME (ng/cm); LB400`, var.equal = TRUE)



exp.mpuf <- exp.data %>%
  filter(sampler == "PUF") %>%
  group_by(time, treatment) %>%
  mutate(PCB4) %>%
  distinct(`PCB4`) %>%
  rename("PCB 4 mass in PUF (ng)" = `PCB4`)




t.test(exp.mpuf$`PCB 4 mass in PUF (ng)`, exp.mspme$`PCB 4 mass in SPME (ng/cm)`)





PCB.4 <- cbind(exp.mspme.control, exp.mspme.LB400$`PCB 4 mass in SPME (ng/cm); LB400`,
               exp.mpuf.control$`PCB 4 mass in PUF (ng); Control`, exp.mpuf.LB400$`PCB 4 mass in PUF (ng); LB400`)
pCB.4 <-   add_row(time = 0, "PCB 4 mass in SPME (ng/cm); Control" = 0,
                   "PCB 4 mass in PUF (ng); Control" = 0, 
                   "PCB 4 mass in SPME (ng/cm); LB400" = 0,
                   "PCB 4 mass in PUF (ng); LB400" = 0, .before = 1)

a <- inner_join(exp.mpuf.control, exp.mpuf.LB400, by = "time")


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

