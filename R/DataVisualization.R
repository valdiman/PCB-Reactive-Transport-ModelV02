# Statistical analysis to review data, mostly between control
# and experiments, same time points.

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

# SPME Calibration Analysis -----------------------------------------------
# Read Data

# Read data from Excel sheets
data.1 <- read_excel("Data/SPMECalibration.xlsx", sheet = "day5_14")
data.2 <- read_excel("Data/SPMECalibration.xlsx", sheet = "day42")
data.3 <- read_excel("Data/SPMECalibration.xlsx", sheet = "day85")
data.4 <- read_excel("Data/SPMECalibration.xlsx", sheet = "nonshaken")

# Select individual congeners from datasets
pcbi <- "PCB4"

# Extract relevant columns from each dataset
d.1.pcbi <- data.1[, c("sample", "treatment", "replicate",
                       "time", "length", pcbi)]
d.2.pcbi <- data.2[, c("sample", "treatment", "replicate",
                       "time", "length", pcbi)]
d.3.pcbi <- data.3[, c("sample", "treatment", "replicate",
                       "time", "length", pcbi)]
d.4.pcbi <- data.4[, c("sample", "treatment", "replicate",
                       "time", "length", pcbi)]

# Combine data frames into one
pcb.i <- rbind(d.1.pcbi, d.2.pcbi, d.3.pcbi, d.4.pcbi)

# Create a new variable to group replicates r.3.1, r.3.2,
# and r.3.3 into one category
pcb.i$replicate_grouped <- ifelse(pcb.i$replicate %in% c("r.3.1",
                                                         "r.3.2", "r.3.3"),
                                  "r.3.x",
                                  ifelse(pcb.i$replicate == "r.3", "r.3",
                                         pcb.i$replicate))

# Create plot
ggplot(pcb.i, aes(x = time, y = !!as.name(pcbi)/length,
                  color = replicate_grouped, shape = treatment)) +
  geom_point(size = ifelse(pcb.i$treatment %in% c("calibration",
                                                  "nonshaken"), 2.5, 1)) +
  scale_color_manual(values = c("r.1" = "blue", "r.2" = "red",
                                "r.3" = "green", "r.3.x" = "purple")) +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold(.(pcbi) ~ "(ng/cm)")),
       color = "Replicates",
       shape = "Treatment") +
  theme_bw() +
  theme(legend.position = "right")


# Read data ---------------------------------------------------------------
exp.data <- read.csv("Data/PCBDataV02.csv")

# Organize SPME data -----------------------------------------------------------
# spme = SPME fiber sampler [ng/cm]
# Pull congener-specific data from the dataset & calculate mean
# values for each sampler-treatment combination at each time point
# (i) Select congener, i
# PCBs 4, 17, 19, 31, 52
i <- "PCB19"
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

