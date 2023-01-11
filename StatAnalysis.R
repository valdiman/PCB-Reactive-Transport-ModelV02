# Statistical analysis to review data, mostly between control
# and experiments, same time points.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("ggplot2")

# Load libraries
library(dplyr) # organize data
library(ggplot2) # plotting

# Read data ---------------------------------------------------------------
exp.data <- read.csv("PCBDataV02.csv")

# Organize SPME data -----------------------------------------------------------
# spme = SPME fiber sampler [ng/cm]
# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
# (i) Select congener, i
# PCBs 4, 17, 19, 31, 52
i <- "PCB52"
{exp.mspme <- exp.data %>%
  mutate(exp.data[all_of(i)]/length) %>%
  filter(sampler == "SPME") %>%
  group_by(time, treatment) %>%
  select(all_of(i)) %>%
  rename(!!paste(i, ".SPME", sep ="") := `i`)
print(exp.mspme)

# Plot
ggplot(exp.mspme, aes_string(x = colnames(exp.mspme[1]),
                             y = colnames(exp.mspme[3]),
                      group = colnames(exp.mspme[2]))) +
  geom_point(aes(shape = treatment, color = treatment)) +
  # x-axis not scale
  scale_x_continuous(labels = c('16', '', '35', '', '75'))}

# (ii) Organize data for t test
# Select time, t
t <- "1"
{exp.mspme.t <- exp.mspme %>%
  filter(time == t)
exp.mspme.t <- data.frame(Column1 = exp.mspme.t[c(1:3), 3], 
           Column2 = exp.mspme.t[c(4:6), 3])
exp.mspme.t <- cbind(c(t, t, t), exp.mspme.t[1],
                       exp.mspme.t[2])
exp.mspme.t <- as.data.frame(apply(exp.mspme.t, 2,
                                     as.numeric))
colnames(exp.mspme.t) <- c("time", paste(i, ".SPME[Control]", sep = ""),
                           paste(i, ".SPME[LB400]", sep = ""))
print(exp.mspme.t)}

# Perform T test, including variance comparison (> or < 3 times)
if (var(exp.mspme.t[2],na.rm=TRUE)/var(exp.mspme.t[3], na.rm=TRUE > 3) ||
    var(exp.mspme.t[3], na.rm=TRUE)/var(exp.mspme.t[2], na.rm=TRUE > 3)){
  t.test(exp.mspme.t[2], exp.mspme.t[3])
} else {
  t.test(exp.mspme.t[2],
         exp.mspme.t[3], var.equal = TRUE)
}

# Organize PUF data -------------------------------------------------------
# puf = PUF sampler [ng]
# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
# (i) Select congener, i
# PCBs 4, 17, 19, 31, 52
{exp.mpuf <- exp.data %>%
  mutate(exp.data[all_of(i)]) %>%
  filter(sampler == "PUF") %>%
  group_by(time, treatment) %>%
  select(all_of(i)) %>%
  rename(!!paste(i, ".PUF", sep ="") := `i`)
print(exp.mpuf)

# Plot
ggplot(exp.mpuf, aes_string(x = colnames(exp.mpuf[1]),
                             y = colnames(exp.mpuf[3]),
                             group = colnames(exp.mpuf[2]))) +
  geom_point(aes(shape = treatment, color = treatment)) +
  # x-axis not scale
  scale_x_continuous(labels = c('16', '', '35', '', '75'))
}

# (ii) Organize data for t test
# Select time, t
{exp.mpuf.t <- exp.mpuf %>%
  filter(time == t)
exp.mpuf.t <- data.frame(Column1 = exp.mpuf.t[c(1:3), 3], 
                          Column2 = exp.mpuf.t[c(4:6), 3])
exp.mpuf.t <- cbind(c(t, t, t), exp.mpuf.t[1],
                     exp.mpuf.t[2])
exp.mpuf.t <- as.data.frame(apply(exp.mpuf.t, 2,
                                   as.numeric))
colnames(exp.mpuf.t) <- c("time", paste(i, ".PUF[Control]", sep = ""),
                           paste(i, ".PUF[LB400]", sep = ""))
print(exp.mpuf.t)}

# Perform T test, including variance comparison (> or < 3 times)
if (var(exp.mpuf.t[2],na.rm=TRUE)/var(exp.mpuf.t[3], na.rm=TRUE > 3) ||
    var(exp.mpuf.t[3], na.rm=TRUE)/var(exp.mpuf.t[2], na.rm=TRUE > 3)){
  t.test(exp.mpuf.t[2], exp.mpuf.t[3])
} else {
  t.test(exp.mpuf.t[2],
         exp.mpuf.t[3], var.equal = TRUE)
}

