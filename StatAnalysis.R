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
print(exp.mspme)

# Plot
ggplot(exp.mspme, aes(x = time, y = PCB4.SPME, group = treatment)) +
  geom_point(aes(shape = treatment, color = treatment))

# (ii) Organize data for t test
# Select time, t
t <- "3"
exp.mspme.t <- exp.mspme %>%
  filter(time == t)
exp.mspme.t <- data.frame(Column1 = exp.mspme.t$`PCB4.SPME`[c(1:3)], 
           Column2 = exp.mspme.t$`PCB4.SPME`[c(4:6)])
exp.mspme.t <- cbind(c(t, t, t), exp.mspme.t$Column1,
                       exp.mspme.t$Column2)
exp.mspme.t <- as.data.frame(apply(exp.mspme.t, 2,
                                     as.numeric))
colnames(exp.mspme.t) <- c("time", "PCB4.SPME[Control]",
                             "PCB4.SPME[LB400]")
print(exp.mspme.t)

# t-test
t.test(exp.mspme.t.1$`PCB4.SPME[Control]`,
       exp.mspme.t.1$`PCB4.SPME[LB400]`, var.equal = TRUE)
#
t.test(exp.mspme.t.1$`PCB4.SPME[Control]`,
       exp.mspme.t.1$`PCB4.SPME[LB400]`)
